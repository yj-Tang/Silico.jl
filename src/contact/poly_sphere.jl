################################################################################
# contact
################################################################################
struct PolySphere{T,D,NP,NC} <: Contact{T,D,NP,NC}
    name::Symbol
    parent_name::Symbol
    child_name::Symbol
    index::NodeIndices
    parent_shape::Shape{T}
    child_shape::Shape{T}
    friction_coefficient::Vector{T}
end

function PolySphere(parent_body::AbstractBody{T,D}, child_body::AbstractBody{T};
        parent_shape_id::Int=1,
        child_shape_id::Int=1,
        name::Symbol=:contact,
        friction_coefficient=0.2) where {T,D}

    parent_name = parent_body.name
    child_name = child_body.name
    parent_shape = deepcopy(parent_body.shapes[parent_shape_id])
    child_shape = deepcopy(child_body.shapes[child_shape_id])

    index = NodeIndices()

    Np = constraint_dimension(parent_shape)
    Nc = constraint_dimension(child_shape)
    return PolySphere{T,D,Np,Nc}(
        name,
        parent_name,
        child_name,
        index,
        parent_shape,
        child_shape,
        [friction_coefficient],
        )
end

primal_dimension(contact::PolySphere{T,D}) where {T,D} = D # c
# because ϕ can be achieved
cone_dimension(contact::PolySphere{T,D,NP}) where {T,D,NP} = 1 + 1 + 2*(D-1) + NP # γ ψ β λp

function unpack_variables(x::Vector, contact::PolySphere{T,D,NP}) where {T,D,NP}
    num_cone = cone_dimension(contact)
    off = 0

    if D == 2
        c = x[off .+ (1:2)]; off += 2

        γ = x[off .+ (1:1)]; off += 1
        ψ = x[off .+ (1:1)]; off += 1
        β = x[off .+ (1:2)]; off += 2
        λp = x[off .+ (1:NP)]; off += NP
    
        sγ = x[off .+ (1:1)]; off += 1
        sψ = x[off .+ (1:1)]; off += 1
        sβ = x[off .+ (1:2)]; off += 2
        sp = x[off .+ (1:NP)]; off += NP
    elseif D == 3
        nβ = (D - 1) * 2

        c = x[off .+ (1:D)]; off += D

        γ = x[off .+ (1:1)]; off += 1
        ψ = x[off .+ (1:1)]; off += 1
        β = x[off .+ (1:nβ)]; off += nβ
        λp = x[off .+ (1:NP)]; off += NP
        # λc = x[off .+ (1:NC)]; off += NC

        sγ = x[off .+ (1:1)]; off += 1
        sψ = x[off .+ (1:1)]; off += 1
        sβ = x[off .+ (1:nβ)]; off += nβ
        sp = x[off .+ (1:NP)]; off += NP
    end
    return c, γ, ψ, β, λp, sγ, sψ, sβ, sp
end

function split_parameters(θ, contact::PolySphere{T}) where T
    friction_coefficient, parent_parameters, child_parameters = unpack_parameters(θ, contact)
    shape_p = contact.parent_shape
    shape_c = contact.child_shape
    Ap, bp, op = unpack_parameters(shape_p, parent_parameters)
    bop = bp + Ap * op
    radc, offc = unpack_parameters(shape_c, child_parameters)
    return friction_coefficient, Ap, bop, radc, offc
end

# min 0.5 * (x - c)' * (x - c)
# s.t. Ax <= b
# function residual!(e, x, θ, contact::PolySphere{T,2,NP},
#         pbody::AbstractBody, cbody::AbstractBody) where {T,NP}

#     # unpack parameters
#     friction_coefficient, Ap, bop, radc, offc =
#         split_parameters(θ[contact.index.parameters], contact)

#     pp2, timestep_p = unpack_pose_timestep(θ[pbody.index.parameters], pbody)
#     pc2, timestep_c = unpack_pose_timestep(θ[cbody.index.parameters], cbody)

#     # unpack variables
#     c, γ, ψ, β, λp, sγ, sψ, sβ, sp = unpack_variables(x[contact.index.variables], contact)
#     vp25 = unpack_variables(x[pbody.index.variables], pbody)
#     vc25 = unpack_variables(x[cbody.index.variables], cbody)
#     pp3 = pp2 + timestep_p[1] * vp25
#     pc3 = pc2 + timestep_c[1] * vc25

#     # contact position in the world frame
#     contact_w = c + pp3[1:2]
#     # contact_p is expressed in pbody's frame
#     contact_p = x_2d_rotation(pp3[3:3])' * (contact_w - pp3[1:2])
#     # contact_c is expressed in cbody's frame
#     contact_c = x_2d_rotation(pc3[3:3])' * (contact_w - pc3[1:2])

#     # signed distance function
#     ϕ = [sqrt((pc3[1:2] + offc - contact_w)' * (pc3[1:2] + offc - contact_w))] - radc

#     # contact normal and tangent in the world frame
#     normal_pw = -x_2d_rotation(pp3[3:3]) * Ap' * λp
#     normal_cw = contact_w - pc3[1:2] - offc
#     R = [0 1; -1 0]
#     tangent_pw = R * normal_pw
#     tangent_cw = R * normal_cw

#     # rotation matrix from contact frame to world frame
#     wRp = [tangent_pw normal_pw] # n points towards the parent body, [t,n,z] forms an oriented vector basis
#     wRc = [tangent_cw normal_cw] # n points towards the parent body, [t,n,z] forms an oriented vector basis

#     # force at the contact point in the contact frame
#     f = [β[1] - β[2]; γ]
#     # force at the contact point in the world frame
#     f_pw = +wRp * f # parent
#     f_cw = -wRc * f # child
#     # torques at the centers of masses in world frame
#     τ_pw = (skew([contact_w - pp3[1:2]; 0]) * [f_pw; 0])[3:3]
#     τ_cw = (skew([contact_w - pc3[1:2]; 0]) * [f_cw; 0])[3:3]
#     # overall wrench on both bodies in world frame
#     # mapping the contact force into the generalized coordinates (at the centers of masses and in the world frame)
#     wrench_p = [f_pw; τ_pw]
#     wrench_c = [f_cw; τ_cw]

#     # tangential velocities at the contact point
#     tanvel_p = vp25[1:2] + (skew([pp3[1:2] - contact_w; 0]) * [zeros(2); vp25[3]])[1:2]
#     tanvel_p = tanvel_p' * tangent_pw
#     tanvel_c = vc25[1:2] + (skew([pc3[1:2] - contact_w; 0]) * [zeros(2); vc25[3]])[1:2]
#     tanvel_c = tanvel_c' * tangent_cw
#     tanvel = tanvel_p - tanvel_c

#     # contact equality
#     optimality = [
#         (contact_w - pc3[1:2] - offc) + x_2d_rotation(pp3[3:3]) * Ap' * λp;
#     ]
#     slackness = [
#         sγ - ϕ;
#         sψ - (friction_coefficient[1] * γ - [sum(β)]);
#         sβ - ([+tanvel; -tanvel] + ψ[1]*ones(2));
#         sp - (- Ap * contact_p + bop);
#     ]

#     # fill the equality vector (residual of the equality constraints)
#     e[contact.index.optimality] .+= optimality
#     e[contact.index.slackness] .+= slackness
#     e[pbody.index.optimality] .-= wrench_p
#     e[cbody.index.optimality] .-= wrench_c
#     return nothing
# end

function residual!(e, x, θ, contact::PolySphere{T,3,NP},
    pbody::AbstractBody, cbody::AbstractBody) where {T,NP}

# unpack parameters
friction_coefficient, Ap, bop, radc, offc =
    split_parameters(θ[contact.index.parameters], contact)

pp2, timestep_p = unpack_pose_timestep(θ[pbody.index.parameters], pbody)
pc2, timestep_c = unpack_pose_timestep(θ[cbody.index.parameters], cbody)
xp2 = pp2[1:3]
qp2 = pp2[4:7]
xc2 = pc2[1:3]
qc2 = pc2[4:7]

# unpack variables
c, γ, ψ, β, λp, sγ, sψ, sβ, sp = unpack_variables(x[contact.index.variables], contact)
vp25, ϕp25 = unpack_variables(x[pbody.index.variables], pbody)
vc25, ϕc25 = unpack_variables(x[cbody.index.variables], cbody)
xp3 = xp2 + timestep_p[1] * vp25
xc3 = xc2 + timestep_c[1] * vc25
qp3 = quaternion_increment(qp2, timestep_p[1] * ϕp25)
qc3 = quaternion_increment(qc2, timestep_c[1] * ϕc25)


# contact position in the world frame
contact_w = c + xp3
# contact_p is expressed in pbody's frame
contact_p = vector_rotate(dagger(qp3), contact_w - xp3)
# contact_c is expressed in cbody's frame
contact_c = vector_rotate(dagger(qc3), contact_w - xc3)

# signed distance function
ϕ = [sqrt((xc3 + offc - contact_w)' * (xc3 + offc - contact_w))] - radc

# contact normal and tangent in the world frame
normal_pw = -vector_rotate(qp3, Ap' * λp)
normal_cw = contact_w - xc3 - offc
tangent_candidate = [1, 0, 0] # need to use a vector orthogonal to the normal at timestep 2
tangent_pw1, tangent_pw2, wRp = tangential_plane(normal_pw, tangent_candidate=tangent_candidate)
tangent_cw1, tangent_cw2, wRc = tangential_plane(normal_cw, tangent_candidate=tangent_candidate)

# force at the contact point in the contact frame
f = [β[1:2] - β[3:4]; γ]
# force at the contact point in the world frame
f_pw = +wRp * f # parent
f_cw = -wRc * f # child
# torques at the centers of masses in world frame
τ_pw = skew(contact_w - xp3) * f_pw
τ_cw = skew(contact_w - xc3) * f_cw
# overall wrench on both bodies in world frame
# mapping the contact force into the generalized coordinates (at the centers of masses and in the world frame)
torsional_velocity_pw = vector_rotate(qp2, 2ϕp25)' * normal_pw * normal_pw
torsional_velocity_cw = vector_rotate(qc2, 2ϕc25)' * normal_cw * normal_cw
torsional_velocity_w = torsional_velocity_pw - torsional_velocity_cw
torsional_drag = -torsional_velocity_w .* friction_coefficient / 20
wrench_p = [f_pw; vector_rotate(dagger(qp2), τ_pw + torsional_drag)] # we should apply Euler equation (rotational part) in the body frame, not the world frame
wrench_c = [f_cw; vector_rotate(dagger(qc2), τ_cw - torsional_drag)] # we should apply Euler equation (rotational part) in the body frame, not the world frame

# tangential velocities at the contact point
tanvel_p = vp25 + skew(xp3 - contact_w) * vector_rotate(qp2, 2ϕp25) # let's assume 2ϕp25 is the angular velocity
tanvel_c = vc25 + skew(xc3 - contact_w) * vector_rotate(qc2, 2ϕc25) # let's assume 2ϕp25 is the angular velocity
tanvel_p = (1 - tanvel_p' * normal_pw) * tanvel_p
tanvel_c = (1 - tanvel_c' * normal_cw) * tanvel_c
tanvel_p = [tanvel_p' * tangent_pw1, tanvel_p' * tangent_pw2]
tanvel_c = [tanvel_c' * tangent_cw1, tanvel_c' * tangent_cw2]
tanvel = tanvel_p - tanvel_c

# contact equality
optimality = [
    (contact_w - xc3 - offc) + vector_rotate(qp3, Ap' * λp);
]

slackness = [
    sγ - ϕ;
    sψ - (friction_coefficient[1] * γ - [sum(β)]);
    sβ - ([+tanvel; -tanvel] + ψ[1]*ones(4));
    sp - (- Ap * contact_p + bop);
]

# fill the equality vector (residual of the equality constraints)
e[contact.index.optimality] .+= optimality
e[contact.index.slackness] .+= slackness
e[pbody.index.optimality] .-= wrench_p
e[cbody.index.optimality] .-= wrench_c
return nothing
end

