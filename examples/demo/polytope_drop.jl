using Plots
using Statistics
using Random

using Silico: Visualizer, set_floor!, set_light!, set_background!, Diagonal, get_3d_polytope_drop, Mehrotra, simulate!, build_mechanism!, set_mechanism!, visualize!, normalize, planar_sliding, sliding_on_slope, get_polytope_insertion, get_polytope_drop

################################################################################
# visualization
################################################################################
vis = Visualizer()
open(vis)
set_floor!(vis)
set_light!(vis)
set_background!(vis)

################################################################################
# define mechanism
################################################################################
timestep = 0.05;
gravity = -9.81;
mass = 1.0;
inertia = 0.2 * ones(1);

mech = get_polytope_drop(;
    timestep=0.05,
    gravity=-9.81,
    mass=1.0,
    inertia=0.2 * ones(1,1),
    friction_coefficient=0.9,
    # method_type=:symbolic,
    method_type=:finite_difference,
    options=Mehrotra.Options(
        verbose=false,
        complementarity_tolerance=1e-3,
        # compressed_search_direction=true,
        compressed_search_direction=false,
        sparse_solver=false,
        warm_start=false,
        )
    );

# solve!(mech.solver)
################################################################################
# test simulation
################################################################################
xp2 = [+0.0,1.5,-0.001]
vp15 = [-0,0,-0.0]
z0 = [xp2; vp15]

u0 = ones(3)*1
H0 = 150

@elapsed storage = simulate!(mech, z0, H0)

################################################################################
# visualization
################################################################################
build_mechanism!(vis, mech)
set_mechanism!(vis, mech, storage, 10)

visualize!(vis, mech, storage, build=false)

# scatter(storage.iterations)
# plot!(hcat(storage.variables...)')


################################################################################
# planning with gradient 
################################################################################

## initializtion 
# dimensions 
n = 6;
m = 3;
T = 31;

# set initial pose
u0 = zeros(3)
z0 = [0.0, 1.0, 0.0, 0.0, 0.0, 0.0]
z1 = zero(z0)

# define Cost
goal = [0.0, 2.0, 0.0, 0.0, 0.0, 0.0]
Q = Diagonal(1.0 * [1.0, 1.0, 1.0, 1.0, 1.0, 1.0])
P = Diagonal(0.1 * [1.0, 1.0, 1.0])
cost = (z, u) -> transpose(z - goal) * Q * (z - goal) + transpose(u) * P * u

# find gradients
∇c_∇z = (z) -> transpose(2*Q*(z - goal))

# ∇z_∇u: calculated by finite difference
dz0 = zeros(length(z0), length(z0))
Silico.dynamics_jacobian_state(dz0, mech, z0, ū[1])
du0 = zeros(length(z0), length(u0))
Silico.dynamics_jacobian_input(du0, mech, z0, ū[1])
# or
∇z_∇z = jacobian(central_fdm(5, 1), z -> Silico.dynamics_FD(z1, mech, z, ū[1]), z0)[1]
∇z_∇u = jacobian(central_fdm(5, 1), u -> Silico.dynamics_FD(z1, mech, z0, u), u0)[1]

# planning
H = 30
z_history = [z0]
α = 0.1
z = z0
dz = zeros(length(z0), length(z0))
for t in 1:H
    global α
    #  using the gradient from finit difference
    # ∇z_∇z = jacobian(central_fdm(5, 1), z -> Silico.dynamics_FD(z1, mech, z, ū[1]), z0)[1]
    # global z = z - transpose(α*∇c_∇z(z)*∇z_∇z)

    #  using the gradient from Silico
    Silico.dynamics_jacobian_state(dz, mech, z, ū[1])
    global z = z - transpose(α*∇c_∇z(z)*dz)
    push!(z_history, z)
end

visualize!(vis, mech, z_history)


################################################################################
# planning with iLQR
################################################################################
using IterativeLQR

## initializtion 
# dimensions 
n = 6;
m = 3;
T = 31;

# set initial pose
u0 = zeros(3)
z0 = [0.0, 1.0, 0.0, 0.0, 0.0, 0.0]

# rollout
dyn = IterativeLQR.Dynamics(
    (y, z, u, w) -> Silico.dynamics(y, mech, z, u),
    (dz, z, u, w) -> Silico.dynamics_jacobian_state(dz, mech, z, u),
    (du, z, u, w) -> Silico.dynamics_jacobian_input(du, mech, z, u),
    n, n, m)
model = [dyn for t = 1:T-1]

ū = [0.1 * randn(m) for t = 1:T-1]
z̄ = IterativeLQR.rollout(model, z0, ū)

# define Cost
goal = [0.0, 2.0, 0.0, 0.0, 0.0, 0.0]
Q = Diagonal(1.0 * [1.0, 1.0, 1.0, 1.0, 1.0, 1.0])
P = Diagonal(0.1 * [1.0, 1.0, 1.0])
cost_ts = [(z, u) -> transpose(z - goal) * Q * (z - goal) + transpose(u) * P * u for t = 1:T-1]
cost_T = (z, u) -> transpose(z - goal) * Q * (z - goal)

cts = [IterativeLQR.Cost(ot, n, m) for ot in ots]
cT = IterativeLQR.Cost(oT, n, 0)
obj = [cts..., cT]

function goal_con(z, u)
    return z - goal
end

goal_con(z_sol[end], u_sol[end])

con_policyt = IterativeLQR.Constraint()
con_policyT = IterativeLQR.Constraint(goal_con, n, 0)
# con_policyT = IterativeLQR.Constraint()

cons = [[con_policyt for t = 1:T-1]..., con_policyT]


## solver
solver = Solver(model, obj, cons)
initialize_controls!(solver, ū) 
initialize_states!(solver, z̄)

## solve
solve!(solver)

## solution
x_sol, u_sol = get_trajectory(solver)

visualize!(vis, mech, x_sol)

goal_con(x_sol[end], zeros(6))

# TODOS
# use finit difference to achieve the gradient?
# maybe try mujoco to test iLQR
# try other iLQR method