######################################################################
# shape
######################################################################
function build_shape!(vis::Visualizer, shape::PolytopeShape{T,Ng,D};
        collider_color=RGBA(0.2, 0.2, 0.2, 0.8),
        ) where {T,Ng,D}

    A = shape.A
    b = shape.b
    o = shape.o
    if D == 2
        build_2d_polytope!(vis, A, b + A * o, color=collider_color)
    else
        build_polytope!(vis, A, b + A * o, color=collider_color)
    end
    return nothing
end

function build_shape!(vis::Visualizer, shape::SphereShape{T,Ng,D};
        collider_color=RGBA(0.2, 0.2, 0.2, 0.8),
        ) where {T,Ng,D}

    δ = 5e-3 * [0, 0, 1]
    offset = [zeros(3-D); shape.position_offset]
    radius = shape.radius[1]
    setobject!(vis,
        HyperSphere(GeometryBasics.Point(offset...), radius),
        MeshPhongMaterial(color=collider_color));
    setobject!(vis[:fake],
        HyperSphere(GeometryBasics.Point((offset .+ δ)...), radius),
        MeshPhongMaterial(color=RGBA(0.8, 0.8, 0.8, 0.8)));
    return nothing
end

function build_shape!(vis::Visualizer, shape::CapsuleShape{T,Ng,D};
        collider_color=RGBA(0.2, 0.2, 0.2, 0.8),
        ) where {T,Ng,D}

    δ = 5e-3 * [0, 0, 1]
    offset = [zeros(3-D); shape.position_offset]
    offset_left = offset + [0; -shape.segment/2; 0]
    offset_right = offset + [0; +shape.segment/2; 0]
    radius = shape.radius[1]
    segment = shape.segment[1]

    setobject!(vis[:cylinder],
        MeshCat.GeometryBasics.Cylinder(Point(offset_left...), Point(offset_right...), radius),
        MeshPhongMaterial(color=collider_color));
    setobject!(vis[:left],
        HyperSphere(GeometryBasics.Point(offset_left...), radius),
        MeshPhongMaterial(color=collider_color));
    setobject!(vis[:right],
        HyperSphere(GeometryBasics.Point(offset_right...), radius),
        MeshPhongMaterial(color=collider_color));
    return nothing
end


function build_shape!(vis::Visualizer, shape::PaddedPolytopeShape110{T,Ng,D};
        collider_color=RGBA(0.2, 0.2, 0.2, 0.8),
        ) where {T,Ng,D}

    radius = shape.radius[1]
    A = shape.A
    b = shape.b
    o = shape.o

    @assert D == 2
    build_2d_polytope!(vis, A, b + A * o, color=collider_color, thickness=2radius)
    h = RobotVisualizer.Polyhedra.hrep(A, b + A * o)
    p = RobotVisualizer.Polyhedra.polyhedron(h)
    v = RobotVisualizer.Polyhedra.vrep(p)
    nv = size(v.V, 1)
    edges = Vector{Vector{Int}}()
    for i = 1:nv
        vi =  v.V[i,:]
        for j = i+1:nv
            vj = v.V[j,:]
            c = 0.5 * h.A * (vi+vj) - h.b
            (maximum(c) > -1e-5) && push!(edges, [i,j])
        end
    end

    material = MeshPhongMaterial(color=collider_color)
    for edge in edges
        vi = v.V[edge[1],:]
        vj = v.V[edge[2],:]
        cylinder = MeshCat.Cylinder(GeometryBasics.Point(0, vi...), GeometryBasics.Point(0, vj...), radius)
        sphere = MeshCat.HyperSphere(GeometryBasics.Point(0, vi...), radius)
        setobject!(vis[Symbol(:cylinder_, edge[1], :_, edge[2])], cylinder, material)
    end
    for i = 1:nv
        vi = v.V[i,:]
        sphere = MeshCat.HyperSphere(GeometryBasics.Point(0, vi...), radius)
        setobject!(vis[Symbol(:sphere_, i)], sphere, material)
    end
    return nothing
end

function build_contact_shape!(vis::Visualizer, shape::Shape; collider_color=nothing)
end

function build_contact_shape!(vis::Visualizer, shape::PolytopeShape;
        collider_color=RGBA(0.5,0.5,0.5,0.0))
    build_shape!(vis, shape; collider_color=collider_color)
    return nothing
end

######################################################################
# body
######################################################################
function build_body!(vis::Visualizer, body::AbstractBody;
        name::Symbol=body.name,
        collider_color=RGBA(0.2, 0.2, 0.2, 0.8),
        center_of_mass_color=RGBA(1, 1, 1, 1.0),
        center_of_mass_radius=0.025,
        ) where T

    # shapes
    for (i,shape) in enumerate(body.shapes)
        build_shape!(vis[:bodies][name][Symbol(i)], shape, collider_color=collider_color)
    end

    # center of mass
    setobject!(vis[:bodies][name][:com],
        HyperSphere(GeometryBasics.Point(0,0,0.), center_of_mass_radius),
        MeshPhongMaterial(color=center_of_mass_color));
    return nothing
end

function set_body!(vis::Visualizer, body::AbstractBody{T,D}, pose; name=body.name) where {T,D}
    x = [0; pose[1:2]]
    q = RotX(pose[3])
    if D == 3
        x = pose[1:3]
        q = Quaternion(pose[4:7]...)
    end
    settransform!(vis[:bodies][name], MeshCat.compose(
        MeshCat.Translation(SVector{3}(x)),
        MeshCat.LinearMap(rotationmatrix(q)),
        )
    )
    return nothing
end

######################################################################
# mechanism
######################################################################
function build_mechanism!(vis::Visualizer, mechanism::Mechanism;
        show_contact::Bool=true,
        color=RGBA(1, 1, 1, 0.8),
        center_of_mass_color=RGBA(1,1,1,1.0),
        env_color=RGBA(0.5, 0.5, 0.5, 0.8),
        name::Symbol=:robot)

    for body in mechanism.bodies
        build_body!(vis[name], body, collider_color=color, center_of_mass_color=center_of_mass_color)
    end
    if show_contact
        for contact in mechanism.contacts
            RobotVisualizer.build_frame!(vis[name][:contacts], dimension=space_dimension(contact), name=contact.name)
        end
    end
    for contact in mechanism.contacts
        build_contact_shape!(vis[contact.name], contact.child_shape; collider_color=env_color)
    end
    return nothing
end

function set_mechanism!(vis::Visualizer, mechanism::Mechanism, storage::TraceStorage,
        i::Int; show_contact::Bool=true, name::Symbol=:robot)

    for (j,body) in enumerate(mechanism.bodies)
        set_body!(vis[name], body, storage.x[i][j])
    end
    if show_contact
        for (j, contact) in enumerate(mechanism.contacts)
            ii = max(1,i-1) # needed otherwise the contact frame is a one step ahead of the bodies.
            origin = storage.contact_point[ii][j]
            normal = storage.normal[ii][j]
            tangent_x = storage.tangent_x[ii][j]
            tangent_y = storage.tangent_y[ii][j]
            RobotVisualizer.set_frame!(vis[name][:contacts], origin, normal, tangent_x, tangent_y, name=contact.name, normalize=true)
        end
    end
    return nothing
end

function set_mechanism!(vis::Visualizer, mechanism::Mechanism, z; name::Symbol=:robot)
    off = 0
    for (j,body) in enumerate(mechanism.bodies)
        nz = state_dimension(body)
        set_body!(vis[name], body, z[off .+ (1:nz)]); off += nz
    end
    return nothing
end

function visualize!(vis::Visualizer, mechanism::Mechanism, storage::TraceStorage{T,H};
        build::Bool=true,
        show_contact::Bool=true,
        color=RGBA(0.2, 0.2, 0.2, 0.8),
        env_color=RGBA(0.5, 0.5, 0.5, 1.0),
        name::Symbol=:robot,
        animation=MeshCat.Animation(Int(floor(1/mechanism.bodies[1].timestep[1])))) where {T,H}

    build && build_mechanism!(vis, mechanism, show_contact=show_contact,
        color=color, env_color=env_color, name=name)
    for i = 1:(1)
        atframe(animation, i) do
            set_mechanism!(vis, mechanism, storage, i, show_contact=show_contact, name=name)
        end
        sleep(1)
    end
    MeshCat.setanimation!(vis, animation)
    return vis, animation
end

function visualize!(vis::Visualizer, mechanism::Mechanism, z;
        build::Bool=true,
        show_contact::Bool=true,
        color=RGBA(0.2, 0.2, 0.2, 0.8),
        env_color=RGBA(0.5, 0.5, 0.5, 1.0),
        name::Symbol=:robot,
        animation=MeshCat.Animation(Int(floor(1/mechanism.bodies[1].timestep[1])))) where {T}

    H = length(z)
    build && build_mechanism!(vis, mechanism, show_contact=show_contact,
        color=color, env_color=env_color, name=name)
    for i = 1:H
        atframe(animation, i) do
            set_mechanism!(vis, mechanism, z[i], name=name)
        end
    end
    MeshCat.setanimation!(vis, animation)
    return vis, animation
end

# function set_frame!(vis::Visualizer, origin, normal, tangent_x, tangent_y;
#         name::Symbol=:contact,
#         normalize::Bool=true,
#         axis_length=0.15)
#     dimension = length(origin)
#     if dimension == 2
#         origin = [0; origin]
#         tangent_x = [1; tangent_x]
#         tangent_y = [0; tangent_y]
#         normal = [0; normal]
#     end
#     if normalize
#         tangent_x = axis_length .* tangent_x ./ (norm(tangent_x) + 1e-10)
#         tangent_y = axis_length .* tangent_y ./ (norm(tangent_y) + 1e-10)
#         normal = axis_length .* normal ./ (norm(normal) + 1e-10)
#     end
#
#     settransform!(vis[name][:origin],
#         MeshCat.Translation(MeshCat.SVector{3}(origin...)))
#     set_segment!(vis[name], origin, origin+tangent_x; name=:tangent_x)
#     set_segment!(vis[name], origin, origin+tangent_y; name=:tangent_y)
#     set_segment!(vis[name], origin, origin+normal; name=:normal)
#     return nothing
# end
#
# function build_frame!(vis::Visualizer;
#     dimension::Int=3,
#     name::Symbol=:contact,
#     origin_color=RGBA(0.2, 0.2, 0.2, 0.8),
#     tangent_x_axis_color=RGBA(1, 0, 0, 0.8),
#     tangent_y_axis_color=RGBA(0, 1, 0, 0.8),
#     normal_axis_color=RGBA(0, 0, 1, 0.8),
#     origin_radius=0.025,
#     ) where T
#
#     # axes
#     if dimension == 3
#         build_segment!(vis[name];
#             color=tangent_x_axis_color,
#             segment_radius=origin_radius/2,
#             name=:tangent_x)
#     end
#     build_segment!(vis[name];
#         color=tangent_y_axis_color,
#         segment_radius=origin_radius/2,
#         name=:tangent_y)
#
#     build_segment!(vis[name];
#         color=normal_axis_color,
#         segment_radius=origin_radius/2,
#         name=:normal)
#
#     # origin
#     setobject!(vis[name][:origin],
#         HyperSphere(GeometryBasics.Point(0,0,0.), origin_radius),
#         MeshPhongMaterial(color=origin_color));
#     return nothing
# end



# function RobotVisualizer.build_2d_polytope!(vis::Visualizer, A::Matrix{T}, b::Vector{T};
#         name::Symbol=:polytope,
#         thickness=0.10,
#         color=RGBA(0.8, 0.8, 0.8, 1.0)) where T
#
#     n = size(A)[1]
#     Ae = [zeros(n) A]
#     Ae = [Ae;
#          -1 0 0;
#           1 0 0]
#     be = [b; thickness/2; thickness/2]
#     build_polytope!(vis, Ae, be, name=name, color=color)
#     return nothing
# end
