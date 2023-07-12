function planar_pushing(;
    timestep=0.05,
    gravity=-9.81,
    mass=1.0,
    inertia=0.2 * ones(1,1),
    friction_coefficient=0.9,
    A=[
        +0 +0 +1;
        +0 +0 -1;
        +0 +1 +0;
        +0 -1 +0;
        +1 +0 +0;
        -1 +0 +0;
        ],
    b=0.25*[1,1,1,1,1,1],
    method_type::Symbol=:finite_difference,
    options=Mehrotra.Options(
        complementarity_tolerance=1e-4,
        compressed_search_direction=false,
        max_iterations=30,
        sparse_solver=false,
        warm_start=true,
        )
    )

    # nodes
    object_shapes = [PolytopeShape(A, b)]
    robot_shapes = [SphereShape(0.2, zeros(3))]
    bodies = [
        Body(timestep, mass, inertia, object_shapes, gravity=+gravity, name=:pbody, D=3),
        GroundRobot(timestep, mass, inertia, robot_shapes, gravity=gravity, stiffness=[1e+2, 1e+2, 1e+2], name=:robot)
        ]
    normal = [0.0, 0.0, 1.0]
    floor_shape = HalfspaceShape(normal)
    # PolySphere only support the 2d case.
    # TODO: add 3d case for the PolySphere function
    contacts = [
        PolyHalfSpace(bodies[1], floor_shape;
            name=:floor,
            friction_coefficient=friction_coefficient),
        PolySphere(bodies[1], bodies[2],
        friction_coefficient=friction_coefficient,
        name=:object_robot),
        SphereHalfSpace(bodies[2], floor_shape,
        friction_coefficient=friction_coefficient,
        name=:robot_floor)
        ]
    indexing!([bodies; contacts])

    local_mechanism_residual(primals, duals, slacks, parameters) =
        mechanism_residual(primals, duals, slacks, parameters, bodies, contacts)

    mechanism = Mechanism(
        local_mechanism_residual,
        bodies,
        contacts,
        options=options,
        method_type=method_type)

    Mehrotra.initialize_solver!(mechanism.solver)
    return mechanism
end
