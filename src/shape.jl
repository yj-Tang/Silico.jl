abstract type Shape{T,Ng,D} end

space_dimension(shape::Shape{T,Ng,D}) where {T,Ng,D} = D
constraint_dimension(shape::Shape{T,Ng}) where {T,Ng} = Ng
function constraint_jacobian_β(shape::Shape{T,Ng}, p, α, β) where {T,Ng}
    return zeros(T,Ng,0)
end


################################################################################
# polytope shape
################################################################################
struct PolytopeShape{T,Ng,D} <: Shape{T,Ng,D}
    A::Matrix{T}
    b::Vector{T}
    o::Vector{T}
end

function PolytopeShape(A, b::Vector{T}, o=zeros(T,size(A,2))) where T
    Ng = length(b)
    D = size(A, 2)
    return PolytopeShape{T,Ng,D}(A, b, o)
end

primal_dimension(shape::PolytopeShape) = 0
cone_dimension(shape::PolytopeShape{T,Ng}) where {T,Ng} = Ng
parameter_dimension(shape::PolytopeShape{T,Ng,D}) where {T,Ng,D} = Ng + D * (Ng + 1)
get_parameters(shape::PolytopeShape) = [vec(shape.A); shape.b; shape.o]

function set_parameters!(shape::PolytopeShape{T,Ng,D}, parameters) where {T,Ng,D}
    off = 0
    shape.A .= reshape(parameters[off .+ (1:D*Ng)], (Ng,D)); off += D*Ng
    shape.b .= parameters[off .+ (1:Ng)]; off += Ng
    shape.o .= parameters[off .+ (1:D)]; off += D
    return nothing
end

function unpack_parameters(shape::PolytopeShape{T,Ng,D}, parameters) where {T,Ng,D}
    off = 0
    A = reshape(parameters[off .+ (1:D*Ng)], (Ng,D)); off += D*Ng
    b = parameters[off .+ (1:Ng)]; off += Ng
    o = parameters[off .+ (1:D)]; off += D
    return A, b, o
end

function constraint(shape::PolytopeShape, p, α, β)
    A = shape.A
    b = shape.b
    o = shape.o
    return - A * (p - o) + α .* b
end

function constraint_jacobian_α(shape::PolytopeShape, p, α, β)
    return [shape.b;;]
end

function constraint_jacobian_p(shape::PolytopeShape, p, α, β)
    A = shape.A
    return - A
end

function constraint_jacobian_o(shape::PolytopeShape, p, α, β)
    A = shape.A
    return A
end

################################################################################
# halfspace shape
################################################################################
struct HalfspaceShape{T,Ng,D} <: Shape{T,Ng,D}
    normal::Vector{T}
    position_offset::Vector{T}
end

function HalfspaceShape(normal::AbstractVector{T}, position_offset=zeros(T,length(normal))) where T
    D = length(normal)
    return HalfspaceShape{T,1,D}(normal, position_offset)
end

primal_dimension(shape::HalfspaceShape) = 0
cone_dimension(shape::HalfspaceShape) = 1
parameter_dimension(shape::HalfspaceShape{T,Ng,D}) where {T,Ng,D} = 2D
get_parameters(shape::HalfspaceShape) = [shape.normal; shape.position_offset]

function set_parameters!(shape::HalfspaceShape{T,Ng,D}, parameters) where {T,Ng,D}
    off = 0
    shape.normal .= parameters[off .+ (1:D)]; off += D
    shape.position_offset .= parameters[off .+ (1:D)]; off += D
    return nothing
end

function unpack_parameters(shape::HalfspaceShape{T,Ng,D}, parameters) where {T,Ng,D}
    off = 0
    normal = parameters[off .+ (1:D)]; off += D
    position_offset = parameters[off .+ (1:D)]; off += D
    return normal, position_offset
end

function constraint(shape::HalfspaceShape, p, α, β)
    n = shape.normal
    o = shape.position_offset
    # there is a minus sign because the feasible space for the contact point is 'below' the halfspace
    return [-n' * (p - o)]
end

function constraint_jacobian_α(shape::HalfspaceShape, p, α, β)
    return [0;;]
end

function constraint_jacobian_p(shape::HalfspaceShape, p, α, β)
    n = shape.normal
    return [-n';;]
end

function constraint_jacobian_o(shape::HalfspaceShape, p, α, β)
    n = shape.normal
    return [n';;]
end

################################################################################
# sphere shape
################################################################################
struct SphereShape{T,Ng,D} <: Shape{T,Ng,D}
    radius::Vector{T}
    position_offset::Vector{T}
end

function SphereShape(radius::T, position_offset=zeros(T,2)) where T
    D = length(position_offset)
    return SphereShape{T,1,D}([radius], position_offset)
end

primal_dimension(shape::SphereShape) = 0
cone_dimension(shape::SphereShape) = 1
parameter_dimension(shape::SphereShape{T,Ng,D}) where {T,Ng,D} = 1 + D
get_parameters(shape::SphereShape) = [shape.radius; shape.position_offset]

function set_parameters!(shape::SphereShape{T,Ng,D}, parameters) where {T,Ng,D}
    off = 0
    shape.radius .= parameters[off .+ (1:1)]; off += 1
    shape.position_offset .= parameters[off .+ (1:D)]; off += D
    return nothing
end

function unpack_parameters(shape::SphereShape{T,Ng,D}, parameters) where {T,Ng,D}
    off = 0
    radius = parameters[off .+ (1:1)]; off += 1
    position_offset = parameters[off .+ (1:D)]; off += D
    return radius, position_offset
end

function constraint(shape::SphereShape{T}, p, α, β) where {T}
    r = shape.radius[1]
    o = shape.position_offset
    return [- (p - o)' * (p - o) + α^2 * r^2]
end

function constraint_jacobian_α(shape::SphereShape, p, α, β)
    r = shape.radius[1]
    return [2 * α * r^2;;]
end

function constraint_jacobian_p(shape::SphereShape, p, α, β)
    o = shape.position_offset
    return [-2 * (p - o)';;]
end

function constraint_jacobian_o(shape::SphereShape, p, α, β)
    o = shape.position_offset
    return [2 * (p - o)';;]
end
