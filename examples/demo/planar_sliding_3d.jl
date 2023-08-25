using Plots
using Statistics
using Random

using Silico: Visualizer, set_floor!, set_light!, set_background!, Diagonal, get_3d_polytope_drop, Mehrotra, simulate!, build_mechanism!, set_mechanism!, visualize!, normalize, planar_sliding, sliding_on_slope, get_polytope_insertion, get_polytope_drop, dynamics_jacobian_state, dynamics_jacobian_input

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


mech = planar_sliding(;
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
qp2 = normalize([1,0,0,0.0])
qp2 = normalize([1.,0,0,0.0])
xp2 =  [+0.00; +0.00; +0.12500; qp2]
vp15 = [+0.00, 0.50, +0.00, +0.0,+0.0,+0.0]
z0 = [xp2; vp15]

H0 = 16

mech.solver.solution.primals[1:6] .= deepcopy(vp15)
mech.solver.solution.primals[7:9] .= zeros(3)

# @elapsed storage = simulate!(mech, deepcopy(z0), H0)

################################################################################
# visualization
################################################################################
# build_mechanism!(vis, mech)
# set_mechanism!(vis, mech, storage, 1)

# visualize!(vis, mech, storage, build=false)


# scatter(mech.solver.data.residual.all)
# scatter(mech.solver.data.residual.primals)
# scatter(mech.solver.data.residual.duals)
# scatter(mech.solver.data.residual.slacks)
# scatter(storage.iterations)

# mech.solver.dimensions.primals


################################################################################
# test gradient 
################################################################################
# ∇z_∇u: calculated by finite difference
m = 6
T = 30
ū = [0.1 * randn(m) for t = 1:T-1]

dz0 = zeros(length(z0), length(z0))
dynamics_jacobian_state(dz0, mech, z0, ū[1])
du0 = zeros(length(z0), length(u0))
dynamics_jacobian_input(du0, mech, z0, ū[1])
# or
∇z_∇z = jacobian(central_fdm(5, 1), z -> Silico.dynamics_FD(z1, mech, z, ū[1]), z0)[1]
∇z_∇u = jacobian(central_fdm(5, 1), u -> Silico.dynamics_FD(z1, mech, z0, u), u0)[1]