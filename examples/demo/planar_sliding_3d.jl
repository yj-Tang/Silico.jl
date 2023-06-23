using Plots
using Statistics
using Random

# include("/environment/planar_sliding.jl")


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
vp15 = [+0.00, 2.0, +0.00, +0.0,+0.0,+0.0]
z0 = [xp2; vp15]

H0 = 16

mech.solver.solution.primals[1:6] .= deepcopy(vp15)
mech.solver.solution.primals[7:9] .= zeros(3)

@elapsed storage = simulate!(mech, deepcopy(z0), H0)

################################################################################
# visualization
################################################################################
build_mechanism!(vis, mech)
set_mechanism!(vis, mech, storage, 1)

visualize!(vis, mech, storage, build=false)


scatter(mech.solver.data.residual.all)
scatter(mech.solver.data.residual.primals)
scatter(mech.solver.data.residual.duals)
scatter(mech.solver.data.residual.slacks)
scatter(storage.iterations)

mech.solver.dimensions.primals

# RobotVisualizer.convert_frames_to_video_and_gif("polytope_drop_more_stable")
