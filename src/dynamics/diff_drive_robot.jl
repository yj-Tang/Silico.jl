################################################################################
# body
################################################################################
# TODO 1: change the [x,z,\theta] pose to [x,y,\theta]. To make it a ground robot.  2D in this setting denotes pose [x,z,\theta], so we need 3D here for the ground robot
# TODO 2: change input of the quasistatic robot from te pose difference to the velocity. Add one more property in the struct "velocity", refer to the body struct
struct GroundRobot{T,D} <: AbstractBody{T,D}
    name::Symbol
    index::NodeIndices
    pose::Vector{T}
    input::Vector{T}
    gravity::Vector{T}
    timestep::Vector{T}
    mass::Vector{T}
    inertia::Matrix{T}
    stiffness::Vector{T}
    shapes::Vector
end

function GroundRobot(timestep::T, mass, inertia::Matrix,
        shapes::Vector;
        gravity=-9.81,
        # stiffness=1/timestep^2*[mass, mass, inertia[1]],
        # stiffness=timestep^5*[mass, mass, inertia[1]],
        # stiffness=timestep^5*[mass, mass, inertia[1]],
        # stiffness=1e-1 * mass * gravity * ones(3),
        stiffness=[1e+2, 1e+2, 1e+2],
        name::Symbol=:body,
        index::NodeIndices=NodeIndices()) where T

    D = 3
    return GroundRobot{T,D}(
        name,
        index,
        zeros(7),  # pose
        zeros(6),  # input: velocity
        [gravity],
        [timestep],
        [mass],
        inertia,
        stiffness,
        shapes,
    )
end

primal_dimension(body::GroundRobot{T,D}) where {T,D} = 6
cone_dimension(body::GroundRobot{T,D}) where {T,D} = 0

function parameter_dimension(body::GroundRobot{T,D}) where {T,D}
    @assert D == 3
    nq = 7 # configuration
    nu = 6 # input
    n_gravity = 1 # mass
    n_timestep = 1 # mass
    n_mass = 1 # mass
    n_inertia = 3 # inertia
    n_stiffness = 3  # To be changed
    nθ = nq + nu + n_gravity + n_timestep + n_mass + n_inertia + n_stiffness
    return nθ
end

function unpack_variables(x::Vector, body::GroundRobot{T}) where T
    v25 = x[1:3]
    ϕ25 = x[4:6]
    return v25, ϕ25
end

function get_parameters(body::GroundRobot{T,D}) where {T,D}
    @assert D == 3
    pose = body.pose
    input = body.input

    gravity = body.gravity
    timestep = body.timestep
    mass = body.mass
    inertia = body.inertia
    inertia_vec = diag(inertia)
    stiffness = body.stiffness
    θ = [pose; input; gravity; timestep; mass; inertia_vec; stiffness]
    return θ
end

function set_parameters!(body::GroundRobot{T,D}, θ) where {T,D}
    pose, input, timestep, gravity, mass, inertia, stiffness = unpack_parameters(θ, body)
    body.pose .= pose
    body.input .= input

    body.gravity .= gravity
    body.timestep .= timestep
    body.mass .= mass
    body.inertia .= inertia
    body.stiffness .= stiffness
    return nothing
end

function unpack_parameters(θ::Vector, body::GroundRobot{T,D}) where {T,D}
    @assert D == 3
    off = 0
    pose = θ[off .+ (1:7)]; off += 7
    input = θ[off .+ (1:6)]; off += 6

    gravity = θ[off .+ (1:1)]; off += 1
    timestep = θ[off .+ (1:1)]; off += 1
    mass = θ[off .+ (1:1)]; off += 1
    inertia = Diagonal(θ[off .+ (1:3)]); off += 3
    stiffness = θ[off .+ (1:3)]; off += 3
    return pose, input, timestep, gravity, mass, inertia, stiffness
end
parameter_state_indices(body::GroundRobot) = Vector(1:3)
parameter_input_indices(body::GroundRobot) = Vector(4:6)

function unpack_pose_timestep(θ::Vector, body::GroundRobot{T,D}) where {T,D}
    pose, input, timestep, gravity, mass, inertia, stiffness = unpack_parameters(θ, body)
    return pose, timestep
end

function residual!(e, x, θ, body::GroundRobot)
    index = body.index
    # variables = primals = velocity
    v25, ϕ25 = unpack_variables(x[index.variables], body)
    # parameters
    p2, u, timestep, gravity, mass, inertia, stiffness = unpack_parameters(θ[index.parameters], body)

    # x2 = p2[1:3]
    # q2 = p2[4:7]

    # dynamics
    linear_optimality = v25 - u[1:3];
    # because the ground robot only moves on the ground, so its quarternion is [cos(theta), 0, 0, sin(theta)]
    angular_optimality = ϕ25 - u[4:6]



    # # integrator
    # p3 = p2 + timestep[1] * v25

    # # mass matrix
    # K = Diagonal(stiffness)
    # # dynamics 
    # # u=\delta pose
    # # pose = [x, z, \theta]
    # optimality = timestep[1] * K * (p3 - u) - timestep[1] * [0; mass .* gravity; 0];
    e[index.optimality] .+= [linear_optimality; angular_optimality]
    return nothing
end

function get_current_state(body::GroundRobot{T}) where T
    nz = length(body.pose)

    off = 0
    z = zeros(T,nz)
    z[off .+ (1:nz)] .= body.pose; off += nz
    return z
end

function set_current_state!(body::GroundRobot, z)
    np = pose_dimension(body)

    off = 0
    body.pose .= z[off .+ (1:np)]; off += np
    return nothing
end

function get_next_state!(z, variables, body::GroundRobot{T}) where T
    x2 = body.pose[1:3]
    q2 = body.pose[4:7]
    timestep = body.timestep
    v25, ϕ25 = unpack_variables(variables[body.index.variables], body)
    # ϕ15 = body.velocity[4:6]
    # Δt = timestep[1]
    # Δϕ15 = Δt * ϕ15
    # Δϕ25 = Δt * ϕ25
    # inertia = body.inertia

    # nv = velocity_dimension(body)
    off = 0
    z[off .+ (1:3)] .= x2 + timestep[1] .* v25; off += 3
    z[off .+ (1:4)] .= quaternion_increment(q2, timestep[1] .* ϕ25); off += 4
    # z[off .+ (1:nv)] .= [v25; ϕ25]; off += nv
    return nothing
end

state_dimension(body::GroundRobot) = pose_dimension(body) + velocity_dimension(body)
pose_dimension(body::GroundRobot) = 7
velocity_dimension(body::GroundRobot) = 0
input_dimension(body::GroundRobot) = 6
