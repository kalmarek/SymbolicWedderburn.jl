abstract type AbstractClassFunction{T} end # <: AbstractVector{T} ??

Base.eltype(::AbstractClassFunction{T}) where {T} = T

function LinearAlgebra.dot(
    χ::AbstractClassFunction{T},
    ψ::AbstractClassFunction{T},
) where {T}

    val = sum(length(cc) * χ[i] * ψ[-i] for (i, cc) in enumerate(conjugacy_classes(χ)))

    orderG = sum(length, conjugacy_classes(χ))
    val = div(val,orderG)
    return val
end

function (χ::AbstractClassFunction)(g::PermutationGroups.AbstractPerm)
    for (i, cc) in enumerate(conjugacy_classes(χ))
        g ∈ cc && return χ[i]
    end
    throw(DomainError(g, "element does not belong to conjugacy classes of χ"))
end

####################################
# Characters

mutable struct Character{T,CCl<:AbstractOrbit} <: AbstractClassFunction{T}
    vals::Vector{T}
    inv_of::Vector{Int}
    cc::Vector{CCl}
end

function Character(
    v::AbstractVector{T},
    ccls::AbstractVector{CCl},
) where {T,CCl<:AbstractOrbit}
    χ = Character{T,CCl}(v, _inv_of(ccls), ccls)
    return χ
end

function _inv_of(cc::AbstractVector{<:AbstractOrbit})
    inv_of = zeros(Int, size(cc))
    for (i, c) in enumerate(cc)
        g = inv(first(c))
        inv_of[i] = something(findfirst(k -> g in k, cc), 0)
    end
    any(iszero, inv_of) &&
        throw(ArgumentError("Could not find the conjugacy class for inverse of $(first(cc[findfirst(iszero, inv_of)]))."))
    return inv_of
end

function normalize!(χ::Character)

    id = one(first(first(conjugacy_classes(χ))))
    k = χ(id)

    # computing the degree of χ:
    n = sqrt(inv(dot(χ, χ) * k^2))
    @debug n χ.vals

    # renormalizing χ
    for i in eachindex(χ.vals)
        χ.vals[i] *= n
    end
    return χ
end

PermutationGroups.conjugacy_classes(χ::Character) = χ.cc

Base.@propagate_inbounds function Base.getindex(χ::Character, i::Integer)
    @boundscheck checkbounds(χ.vals, abs(i))
    if i < 0
        return @inbounds χ.vals[χ.inv_of[abs(i)]]
    else
        return @inbounds χ.vals[i]
    end
end
