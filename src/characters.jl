abstract type AbstractClassFunction{T} end # <: AbstractVector{T} ??

Base.eltype(::AbstractClassFunction{T}) where {T} = T

function LinearAlgebra.dot(
    χ::AbstractClassFunction{T},
    ψ::AbstractClassFunction{T},
) where {T}
    # TODO: @assert v.cc == w.cc

    R = parent(χ[1]) # TODO make something better here
    val = zero(R)

    for (i, cc) in enumerate(classes(χ))
        val += R(length(cc)) * χ[i] * ψ[-i]
    end

    orderG = R(sum(length, classes(χ)))
    val *= inv(orderG)
    return val
end

function (χ::AbstractClassFunction)(g::GroupElem)
    for (i, cc) in enumerate(classes(χ))
        g ∈ cc && return χ[i]
    end
    throw(DomainError(g, "element does not belong to conjugacy classes of χ"))
end

####################################
# Characters

mutable struct Character{T,CCl} <: AbstractClassFunction{T}
    vals::Vector{T}
    inv_of::Vector{Int}
    cc::Vector{CCl}

    function Character(
        v::AbstractVector{T},
        ccls::AbstractVector{CCl},
        inv_of = _inv_of(ccls),
    ) where {T,CCl}

        χ = new{T,CCl}(v, inv_of, ccls)

        # initial normalization
        R = parent(first(v))
        id = one(first(first(ccls)))
        if !isone(χ(id))
            χ.vals .*= inv(χ(id))
        end

        # computing the degree of χ:
        deg_χ = sqrt(inv(dot(χ, χ)))
        @debug χ.vals, deg_χ

        # renormalizing χ
        χ.vals .*= deg_χ

        return χ
    end
end

function _inv_of(cc::AbstractVector)
    inv_of = zeros(Int, size(cc))
    for (i, c) in enumerate(cc)
        g = inv(first(c))
        inv_of[i] = something(findfirst(k -> g in k, cc), 0)
    end
    @assert !any(iszero, inv_of) "Could not find the conjugacy class of $g."
    return inv_of
end

classes(χ::Character) = χ.cc

Base.@propagate_inbounds function Base.getindex(χ::Character, i::Integer)
    @boundscheck checkbounds(χ.vals, abs(i))
    if i < 0
        return @inbounds χ.vals[χ.inv_of[abs(i)]]
    else
        return @inbounds χ.vals[i]
    end
end
