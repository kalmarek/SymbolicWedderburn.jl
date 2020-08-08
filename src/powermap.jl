struct PowerMap{T} <: AbstractMatrix{Int}
    cc::Vector{T}
    cache::Matrix{Int}

    function PowerMap(
        ccG::AbstractVector{T},
        exp = exponent(ccG),
    ) where {T<:AbstractOrbit}
        return new{T}(ccG, zeros(Int, length(ccG), exp))
    end
end

Base.axes(pm::PowerMap) = (Base.OneTo(length(pm.cc)), 0:size(pm.cache, 2)-1)
Base.size(pm::PowerMap) = size(pm.cache)

function Base.getindex(p::PowerMap, i::Integer, j::Integer)
    # @boundscheck
    if iszero(p.cache[i, j+1])
        gj = first(p.cc[i])^j
        p.cache[i, j+1] = findfirst(cc -> gj ∈ cc, p.cc)
    end
    return p.cache[i, j+1]
end
