# using Pkg: Pkg
# Pkg.activate("$(@__DIR__)")

using Silico: Visualizer, set_floor!, set_light!, set_background!, Diagonal, get_3d_polytope_drop, Mehrotra, simulate!, build_mechanism!, set_mechanism!, visualize!, normalize, planar_sliding, sliding_on_slope

using Silico

using Plots
using Statistics
using Random

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
gravity = -1*9.81;
mass = 0.2;
inertia = 2.8 * Matrix(Diagonal(ones(3)));
friction_coefficient = 0.5

# A=[
#     +0 +0 +1;
#     +0 +0 -1;
#     +0 +1 +0;
#     +0 -1 +0;
#     +1 +0 +0;
#     -1 +0 +0;
#     +1 +1 +1;
#     +1 +1 -1;
#     +1 -1 +1;
#     +1 -1 -1;
#     -1 +1 +1;
#     -1 +1 -1;
#     -1 -1 +1;
#     -1 -1 -1;
#     ]
# b=0.45*[ones(6); 1.5ones(8)]  # H-representation 
# b=0.45*[ones(6);]

A=[
    +0 +0 +2;
    +0 +0 -2;
    +0 +2 +0;
    +0 -2 +0;
    +1 +0 +0;
    -1 +0 +0;
    ]
b=0.25*[1,1,1,1,1,1]


mech = Silico.planar_pushing(;
    timestep=timestep,
    gravity=gravity,
    mass=mass,
    inertia=inertia,
    friction_coefficient=friction_coefficient,
    A=A,
    b=b,
    # method_type=:symbolic,
    method_type=:finite_difference,
    options=Mehrotra.Options(
        verbose=true,
        complementarity_tolerance=1e-4,
        residual_tolerance=1e-5,
        # compressed_search_direction=true,
        compressed_search_direction=false,
        sparse_solver=false,
        warm_start=true,
        complementarity_backstep=1e-2,
        )
    )

################################################################################
# test simulation
################################################################################

q_object0 = normalize([1,0,0,0.0])
q_robot0 = normalize([1.,0,0,0.0])
x_object0  = [0.00; +0.00; +0.1250;q_object0]
x_robot0 = [-0.70; +0.20; 0.2; q_robot0]
v15_object = [0.0,0.0,0.0,0.0,0.0,0.0]
v15_robot = [0.0,0.0,0.0,0.0,0.0,0.0]

z0 = [x_object0; v15_object; x_robot0]


H0 = 140
# u0 = zeros(9)
x_object_goal  = [0.00, 0.00, 0.00]
x_robot_goal = [-0.30, +1.40, -0.00]

# U = []

# uu=0.
# diata_u=0.001
# for i = 1:H0
#     global uu
#     uu += diata_u

#     α = i/H0
#     u_object = [0,0,0,0,0,0]
#     # u_robot = [uu,0.,0,0,0,0]
#     u_robot = [0.1,0.,0,0,0,0]
#     u = [
#         u_object;
#         u_robot;
#         ]
#     push!(U, u)
# end
# u0 = zeros(12)
# # u0 = [0;0;0; +x_finger10; x_finger20]
# ctrl = Silico.open_loop_controller([u0])
# ctrl = Silico.open_loop_controller(U)

# @elapsed storage = Silico.simulate!(mech, z0, H0, controller=ctrl)

# ################################################################################
# # visualization
# ################################################################################
# build_mechanism!(vis, mech)
# set_mechanism!(vis, mech, storage, 1)

# visualize!(vis, mech, storage, build=false)


# scatter(mech.solver.data.residual.all)
# scatter(mech.solver.data.residual.primals)
# scatter(mech.solver.data.residual.duals)
# scatter(mech.solver.data.residual.slacks)
# scatter(storage.iterations)

# mech.solver.dimensions.primals

# RobotVisualizer.convert_frames_to_video_and_gif("polytope_drop_more_stable")


################################################################################
# test the gradient
################################################################################






################################################################################
# planning with iLQR
################################################################################
using IterativeLQR

# ## dimensions 
n = 20; # number of system states
m = 12; # number of control inputs
T = 31; # horizon

# TODO: need to update Silico.dynamics_jacobian_state/dynamics_jacobian_input to 3D case
dyn = IterativeLQR.Dynamics(
    (y, z, u, w) -> Silico.dynamics(y, mech, z, u),
    (dz, z, u, w) -> Silico.dynamics_jacobian_state(dz, mech, z, u),
    (du, z, u, w) -> Silico.dynamics_jacobian_input(du, mech, z, u),
    n, n, m)  # number of next state, number of state, number of control


model = [dyn for t = 1:T-1]


## rollout
Ū = []
for i = 1:T-1
    u_object = [0.,0,0,0,0,0]
    u_robot = [0.,0.,0,0,0,0]
    u = [
        u_object;
        u_robot;
        ]
    push!(Ū, u)
end

q_object0 = normalize([1,0,0,0.0])
q_robot0 = normalize([1.,0,0,0.0])
x_object0  = [0.00; +0.00; +0.1250;q_object0]
x_robot0 = [-0.70; +0.20; 0.2; q_robot0]
v15_object = [0.0,0.0,0.0,0.0,0.0,0.0]
v15_robot = [0.0,0.0,0.0,0.0,0.0,0.0]

z0 = [x_object0; v15_object; x_robot0]

z̄ = IterativeLQR.rollout(model, z0, Ū)
visualize!(vis, mech, z̄)

# ## objective
############################################################################
x_object0  = [0.00; +0.00; +0.1250;q_object0]
goal = [0.; 0.5; 0.0]
ots = [(z, u) -> transpose(z[1:3] - goal) * Diagonal([0.0, 0.0, 1.0e-1]) * (z[1:3] - goal) + transpose(u) * Diagonal([1.0e-3, 1.0e-3, 1.0e-3, 1.0e-3, 1.0e-3, 1.0e-3, 1.0e-3, 1.0e-3, 1.0e-3, 1.0e-3, 1.0e-3, 1.0e-3]) * u for t = 1:T-1]
oT = (z, u) -> transpose(z[1:3] - goal) * Diagonal(0.0 * [1.0, 1.0, 1.0]) * (z[1:3] - goal)


cts = [IterativeLQR.Cost(ot, n, m) for ot in ots]
cT = IterativeLQR.Cost(oT, n, 0)
obj = [cts..., cT]

function goal_con(z, u)
    return z[1:3] - goal
end

# goal_con(z_sol[end], u_sol[end])

con_policyt = IterativeLQR.Constraint()
con_policyT = IterativeLQR.Constraint(goal_con, n, 0)
# con_policyT = IterativeLQR.Constraint()

cons = [[con_policyt for t = 1:T-1]..., con_policyT]

# ## solver
options = IterativeLQR.Options(
        line_search=:armijo,
        max_iterations=50,
        max_dual_updates=8,
        # min_step_size=1e-2,
        objective_tolerance=1e-3,
        lagrangian_gradient_tolerance=1e-3,
        constraint_tolerance=1e-3,
        # initial_constraint_penalty=1e-1,
        scaling_penalty=10.0,
        max_penalty=1e6,
        verbose=true)

s = IterativeLQR.Solver(model, obj, cons, options=options)

IterativeLQR.initialize_controls!(s, Ū)
IterativeLQR.initialize_states!(s, z̄)


# ## solve
# local_callback!(solver::IterativeLQR.Solver) = continuation_callback!(solver, mech, visualize=true)
Silico.reset!(mech, residual_tolerance=1e-6, complementarity_tolerance=1e-3)
# TODO: replace the constrained_ilqr_solve! here
# IterativeLQR.solve!(solver)

@time IterativeLQR.constrained_ilqr_solve!(s)

# ## solution
z_sol, u_sol = IterativeLQR.get_trajectory(s)

z_sol[end]
# ## visualize
z_view = [[z_sol[1] for t = 1:15]..., z_sol..., [z_sol[end] for t = 1:15]...]

visualize!(vis, mech, z_view)
