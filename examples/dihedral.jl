import GroupsCore #GroupsCore API

struct DihedralGroup <: GroupsCore.Group
    n::Int
end

struct DihedralElement <: GroupsCore.GroupElement
    n::Int
    reflection::Bool
    id::Int
end

Base.one(G::DihedralGroup) = DihedralElement(G.n, false, 0)

Base.eltype(::DihedralGroup) = DihedralElement
function Base.iterate(G::DihedralGroup, prev::DihedralElement=DihedralElement(G.n, false, -1))
    if prev.id + 1 >= G.n
        if prev.reflection
            return nothing
        else
            next = DihedralElement(G.n, true, 0)
        end
    else
        next = DihedralElement(G.n, prev.reflection, prev.id + 1)
    end
    return next, next
end
Base.IteratorSize(::Type{DihedralGroup}) = Base.HasLength()

GroupsCore.order(::Type{T}, G::DihedralGroup) where {T} = convert(T, 2G.n)
GroupsCore.gens(G::DihedralGroup) =
    [DihedralElement(G.n, false, 1), DihedralElement(G.n, true, 0)]

# Base.rand not needed for our purposes here

Base.parent(g::DihedralElement) = DihedralGroup(g.n)
function Base.:(==)(g::DihedralElement, h::DihedralElement)
    return g.n == h.n && g.reflection == h.reflection && g.id == h.id
end

function Base.inv(el::DihedralElement)
    (el.reflection || iszero(el.id)) && return el
    return DihedralElement(el.n, false, el.n - el.id)
end
function Base.:*(a::DihedralElement, b::DihedralElement)
    a.n == b.n || error("Cannot multiply elements from different Dihedral groups")
    id = mod(a.reflection ? a.id - b.id : a.id + b.id, a.n)
    return DihedralElement(a.n, a.reflection != b.reflection, id)
end

Base.copy(a::DihedralElement) = DihedralElement(a.n, a.reflection, a.id)

# optional functions:
function GroupsCore.order(T::Type, el::DihedralElement)
    el.reflection && return T(2)
    iszero(el.id )&& return T(1)
    return T(div(el.n, gcd(el.n, el.id)))
end
