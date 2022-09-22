include(joinpath(module_dir(), "Pygmalion/Pygmalion.jl"))

vis = Visualizer()
open(vis)
# render(vis)
set_background!(vis)
# set_light!(vis, direction="Negative")
set_light!(vis, direction="Positive")
set_floor!(vis, x = 0.1, color=RGBA(0.4,0.4,0.4,0.4))
iterate_color = RGBA(1,1,0,0.6);

################################################################################
# inertial parameters
################################################################################
timestep = 0.05;
gravity = -9.81;
mass = 1.0;
inertia = 0.2 * ones(1,1);

mech = get_bundle_drop(;
    timestep=timestep,
    gravity=gravity,
    mass=mass,
    inertia=inertia,
    friction_coefficient=0.1,
    method_type=:symbolic,
	A=A_recipient,
	b=bo_recipient,
    options=Mehrotra.Options(
        verbose=false,
        complementarity_tolerance=1e-3,
        compressed_search_direction=true,
        max_iterations=30,
        sparse_solver=false,
        differentiate=false,
        warm_start=false,
        complementarity_correction=0.5,
        )
    );
Mehrotra.solve!(mech.solver)

################################################################################
# test simulation
################################################################################
vp15 = [-0,0,-0.6*9.0]
xp2 = [+0.0,2.00,+0.00]
z0 = [xp2; vp15]

horizon = 25
storage = simulate!(mech, z0, horizon)
vis, anim = visualize!(vis, mech, storage, show_contact=false)

################################################################################
# camera parameters
################################################################################
nβ = 40
eye_positions = [[0.0, 3.0] for i=1:horizon]
angles = [-π/2 .+ Vector(range(+0.2π, -0.2π, length=nβ)) for i = 1:horizon]

poses = [x[1] for x in storage.x]
x = [p[1:2] for p in poses] # position
θ = [p[3] for p in poses] # orientation
bRw = [[cos(θ[i]) sin(θ[i]); -sin(θ[i]) cos(θ[i])] for i=1:horizon]
noise = [1*[0.1, 0.1, 0.25] .* (rand(3) .- 0.5) for x in storage.x]
noisy_poses = [storage.x[i][1] + noise[i]  for i = 1:horizon]


α, αmax, v, e, hit_indices = vectorized_ray(eye_positions, angles,
	A_recipient, b_recipient, o_recipient, poses;
	altitude_threshold=0.01,
	max_length=50.00,
	)
v_noisy, e_noisy = noise_transform(v, e, noise)

d0b = e .+ α' .* v
d0b_noisy = e_noisy .+ α' .* v_noisy
d0b_split = [d0b[:,(i-1)*nβ .+ (1:nβ)] for i=1:horizon]
d0b_noisy_split = [d0b_noisy[:,(i-1)*nβ .+ (1:nβ)] for i=1:horizon]
d0w_split = [x[i] .+ bRw[i]' * d0b_split[i] for i = 1:horizon]

for i = 1:horizon
	build_point_cloud!(vis[:point_cloud], nβ; color=RGBA(0.9,0.1,0.1,1), name=Symbol(i))
end
set_2d_point_cloud!(vis, eye_positions, d0b_noisy_split; name=:point_cloud)
set_2d_point_cloud!(vis, eye_positions, d0b_split; name=:point_cloud)
set_2d_point_cloud!(vis, eye_positions, d0w_split; name=:point_cloud)


################################################################################
# Initialization
################################################################################
nh = 5
polytope_dimensions = [nh,nh,nh]
np = length(polytope_dimensions)

θinit, kmres = parameter_initialization(d0b_noisy[:,hit_indices], polytope_dimensions)
Ainit, binit, oinit = unpack_halfspaces(deepcopy(θinit), polytope_dimensions)
visualize_kmeans!(vis, θinit, polytope_dimensions, d0b_noisy[:,hit_indices], kmres)
setvisible!(vis[:cluster], true)
setvisible!(vis[:initial], true)
setvisible!(vis[:cluster], false)
setvisible!(vis[:initial], false)

################################################################################
# optimization
################################################################################
function pack_poses_halfspace_variables(θ, denoise)
	return [θ; vcat(denoise...)]
end

function unpack_poses_halfspace_variables(vars, polytope_dimensions)
	np = length(polytope_dimensions)
	nθ = 3 * sum(polytope_dimensions) + 2 * np
	nv = length(vars)
	H = Int((nv - nθ) / 3)

	θ = vars[1:nθ]
	denoise = [vars[nθ + 3*(i-1) .+ (1:3)] for i = 1:H]
	return θ, denoise
end

# projection
function local_projection(vars)
	θ, denoise = unpack_poses_halfspace_variables(vars, polytope_dimensions)

	θπ = projection(θ, polytope_dimensions,
		Alims=[-1.00, +1.00],
		blims=[+0.05, +0.60],
		olims=[-3.00, +3.00],
		)
	return pack_poses_halfspace_variables(θ, denoise)
end

denoise_init = [zeros(3) for i = 1:horizon]
vars_init = pack_poses_halfspace_variables(θinit, denoise_init)
unpack_poses_halfspace_variables(vars_init, polytope_dimensions)
local_projection(vars_init)

# loss and gradients
function local_loss(vars)
	θ, denoise = unpack_poses_halfspace_variables(vars, polytope_dimensions)
	A, b, bo = preprocess_halfspaces(θ, polytope_dimensions)

	l = noisy_shape_loss(
		α, αmax, v, e, hit_indices,
		A, b, bo,
		noise, denoise,
		ShapeLossOptions1200(),
		)
	return l
end

function local_grad(vars)
	θ, denoise = unpack_poses_halfspace_variables(vars, polytope_dimensions)
	A, b, bo = preprocess_halfspaces(θ, polytope_dimensions)

	grads = noisy_shape_gradient(
		α, αmax, v, e, hit_indices,
		A, b, bo,
		noise, denoise,
		ShapeLossOptions1200(),
		)

	dldAb = [vcat(vec.(grads[1])...); vcat(grads[2]...); vcat(grads[3]...)]
	dlddenoise = vcat(grads[4]...)
	dAbdθ = ForwardDiff.jacobian(θ -> preprocess_halfspaces(θ, polytope_dimensions, vectorize=true), θ)
	return [dAbdθ' * dldAb; dlddenoise]
end


local_loss(vars_init)
local_loss(vars_sol0)
local_grad(vars_init)



################################################################################
# solve
################################################################################
adam_opt = Adam(vars_init, local_loss, local_grad)
# adam_opt.a = 5e-3
adam_opt.a = 12e-3
max_iterations = 400
visual_iterations = 25

subset = 1:Int(ceil(max_iterations/visual_iterations)):max_iterations
@elapsed vars_sol0, vars_iter0 = adam_solve!(adam_opt, projection=local_projection, max_iterations=max_iterations)
θsol0, denoise_sol0 = unpack_poses_halfspace_variables(vars_sol0, polytope_dimensions)

vis, anim = visualize_iterates!(vis, vars_iter0[subset], polytope_dimensions, eye_positions[1],
 	angles, 1e-4, max_iterations=visual_iterations, color=iterate_color)

plot(hcat(noise...)'[:,1], color=:red)
plot!(hcat(denoise_sol0...)'[:,1], color=:blue)

plot(hcat(noise...)'[:,2], color=:red)
plot!(hcat(denoise_sol0...)'[:,2], color=:blue)

plot(hcat(noise...)'[:,3], color=:red)
plot!(hcat(denoise_sol0...)'[:,3], color=:blue)


################################################################################
# visualize denoising
################################################################################
num_hit = sum(hit_indices)
build_point_cloud!(vis[:hit_cloud], num_hit; color=RGBA(0.1,0.1,0.1,1), name=Symbol(1))
settransform!(vis[:hit_cloud], MeshCat.Translation(0.3, 0.0, 0.0))
for (j,i) in enumerate(subset)
	_, denoise_i = unpack_poses_halfspace_variables(vars_iter0[i], polytope_dimensions)
	vi, ei = noise_transform(v, e, noise - denoise_i)
	di = ei .+ α' .* vi
	di_hit = di[:,hit_indices]
	atframe(anim, j) do
		set_2d_point_cloud!(vis, eye_positions[1:1], [di_hit], name=:hit_cloud)
	end
end
MeshCat.setanimation!(vis, anim)
# convert_frames_to_video_and_gif("u_shape_denoising")

Asol, bsol, osol = unpack_halfspaces(θsol0, polytope_dimensions)




















# ## visualizer
vis = Visualizer()
# render(vis)
open(vis)
RobotVisualizer.set_background!(vis)
RobotVisualizer.set_light!(vis)
RobotVisualizer.set_floor!(vis, x=0.1)
RobotVisualizer.set_camera!(vis, zoom=3.0)

################################################################################
# inertial parameters
################################################################################
mech = get_polytope_drop(;
    timestep=timestep,
    gravity=gravity,
    mass=mass,
    inertia=inertia,
    friction_coefficient=0.1,
	method_type=:symbolic,
	A=deepcopy(A_box[1]),
	b=deepcopy(bo_box[1]),
    options=Mehrotra.Options(
        verbose=false,
        complementarity_tolerance=1e-3,
        compressed_search_direction=true,
        max_iterations=30,
        sparse_solver=false,
        differentiate=false,
        warm_start=false,
        complementarity_correction=0.5,
        )
    );
Mehrotra.solve!(mech.solver)

################################################################################
# test simulation
################################################################################
vp15 = [+3.0, +0.0, +2.0]
xp2  = [+0.0, +0.75, +0.0]
z0 = [xp2; vp15]

storage = simulate!(mech, z0, horizon)
vis, anim = visualize!(vis, mech, storage)

configuration_ref = [[storage.x[i][1]; storage.v[i][1]] for i = 1:horizon]

################################################################################
# initial guess
################################################################################
Ap1 = [
	+1.0 +0.0;
	+0.0 +1.0;
	-1.0 +0.0;
	+0.0 -1.0;
	]
bp1 = 0.75*[
	+1,
	+1,
	+1,
	+1,
	]
friction_coefficient1 = 0.5

ϵ_init = 1.0
x1 = [
    0.0; 1.0; 0.0;
    0.0; 0.0; 0.0;
	0.0; 0.0;
	ϵ_init ;
	ϵ_init ;
	ϵ_init ;
	ϵ_init ; ϵ_init ;
	ϵ_init * ones(nh);
	ϵ_init ;

	0.0; 0.0;
	ϵ_init ;
	ϵ_init ;
	ϵ_init ;
	ϵ_init ; ϵ_init ;
	ϵ_init * ones(nh);
	ϵ_init ;

	0.0; 0.0;
    ϵ_init ;
    ϵ_init ;
    ϵ_init ;
    ϵ_init ; ϵ_init ;
    ϵ_init * ones(nh);
    ϵ_init ;

	friction_coefficient1;
	vec(Asol[1]);
	bsol[1] + Asol[1] * osol[1];
	vec(Asol[2]);
	bsol[2] + Asol[2] * osol[2];
    vec(Asol[3]);
    bsol[3] + Asol[3] * osol[3];
    ]
configuration_init = [[deepcopy(configuration_ref[i][1:6]); x1[7:end]] for i = 1:horizon]





# ## dimensions
nh
# (2 + 1 + 1 + 1 + 2 + nh + 1)
nx = 6 + np * (2 + 1 + 1 + 1 + 2 + nh + 1) + 1 + np*nh*3
nu = 0

num_states = [nx for t = 1:horizon]
num_actions = [nu for t = 1:horizon-1]

mech.bodies[1].index
mech.contacts[1].index
mech.contacts[1]

function contact_constraints_equality_t(x, u, mechanism;
        Ac=[0 1.0],
        bc=[0.0],
        )

	# unpack
	body_variables_3, contact_variables_3, contact_parameters_3 = unpack_contact_implicit_variables(x, mechanism)
	c, ϕ, γ, ψ, β, λp, λc = contact_variables_3[1]
    return [γ; ψ; β; λp; λc] .* poly_halfspace_slackness(x, mechanism, Ac=Ac, bc=bc)
end

function contact_constraints_inequality_t(x, u, mechanism;
        Ac=[0 1.0],
        bc=[0.0],
        )

	# unpack
	body_variables_3, contact_variables_3, contact_parameters_3 = unpack_contact_implicit_variables(x, mechanism)
	c, ϕ, γ, ψ, β, λp, λc = contact_variables_3[1]
	friction_coefficient_3, A3, b3 = contact_parameters_3[1]

	p3, v25, c, ϕ, γ, ψ, β, λp, λc, A3, b3, friction_coefficient3 = unpack_state(x)
    return [
        poly_halfspace_slackness(x, mechanism, Ac=Ac, bc=bc);
        γ;
        ψ;
        β;
        λp;
        λc;
    ]
end


unpack_contact_implicit_variables(rand(nx), mech)
polytope_dynamics(rand(nx), rand(nx), rand(nu), mech)
poly_halfspace_slackness(rand(nx), mech)
contact_constraints_equality_t(rand(nx), rand(nu), mech)
contact_constraints_inequality_t(rand(nx), rand(nu), mech)
################################################################################
# problem data
################################################################################
# ## dynamics_model
dynamics_model = [(y, x, u) -> polytope_dynamics(y, x, u, mech; timestep=timestep) for t = 1:horizon-1]

# ## objective
function obj_t(x, u, t, mechanism)
    J = 0.0
    Δp = 1.00 .* (x[1:3] - configuration_ref[t][1:3])
    Δv = 1.00 .* (x[4:6] - configuration_ref[t][4:6])
    # Δθ = 0.50 .* (x[6 + ] - [vec(Ap1); bp1; friction_coefficient1])
    J += 0.5 * dot(Δp, Δp)
    J += 0.5 * dot(Δv, Δv)
    # J += 0.5 * dot(Δθ, Δθ)
    return J
end

objective = [
    (x,u) -> obj_t(x, u, t, mech) for t = 1:horizon
    ]

# ## constraints
equality_t(x, u) = contact_constraints_equality_t(x, u, mech)
equality = [equality_t for t = 1:horizon]

inequality_t(x, u) = contact_constraints_inequality_t(x, u, mech)
nonnegative = [inequality_t for t = 1:horizon]

# ## solver
solver = Solver(objective, dynamics_model, num_states, num_actions,
    equality=equality,
    nonnegative=nonnegative,
    options=Options()
    );

# ## initialize
state_guess = deepcopy(configuration_init)
action_guess = [zeros(0) for t = 1:horizon-1] # may need to run more than once to get good trajectory
initialize_states!(solver, state_guess)
initialize_actions!(solver, action_guess)

# ## solve
solve!(solver)

# ## solution
x_sol, u_sol = get_trajectory(solver)
p_sol = [x[1:3] for x in x_sol]
p_truth = [x[1:3] for x in configuration_ref]

# ## visualize
friction_coefficientp = mean([x[end-45+1] for x in x_sol])
Ap = mean([x[end-3nh-1+1:end-nh-1] for x in x_sol])
Ap = reshape(Ap, (nh,2))
bp = mean([x[end-nh-1+1:end-1] for x in x_sol])

build_2d_polytope!(vis[:polytope_init], Ap1, bp1, name=:polytope, color=MeshCat.RGBA(1,1,1,0.4))
build_2d_polytope!(vis[:polytope_sol], Ap, bp, name=:polytope, color=MeshCat.RGBA(1,0,0,0.4))
build_2d_polytope!(vis[:polytope_truth], Ap0, bp0, name=:polytope, color=MeshCat.RGBA(0,0,0,1.0))
settransform!(vis[:polytope_init][:polytope], MeshCat.Translation(0.00,0,0))
settransform!(vis[:polytope_sol][:polytope], MeshCat.Translation(0.05,0,0))
settransform!(vis[:polytope_truth][:polytope], MeshCat.Translation(0.10,0,0))
anim = MeshCat.Animation()
for ii = 1:horizon+20
    atframe(anim, ii) do
		i = clamp(ii, 1, horizon)
		set_2d_polytope!(vis, p_truth[i][1:2], p_truth[i][3:3], name=:polytope_init)
		set_2d_polytope!(vis, p_sol[i][1:2], p_sol[i][3:3], name=:polytope_sol)
        set_2d_polytope!(vis, p_truth[i][1:2], p_truth[i][3:3], name=:polytope_truth)
    end
end
setanimation!(vis, anim)
u_sol
x_sol
