abstract type AbstractActionHomomorphism{T} end

struct InducingHomomorphism{T,V} <: AbstractActionHomomorphism{T}
    features::V
    reversef::Dict{T,Int}
end

Base.getindex(ihom::InducingHomomorphism, i::Integer) = ihom.features[i]
Base.getindex(ihom::InducingHomomorphism{T}, f::T) where {T} = ihom.reversef[f]
(ihom::InducingHomomorphism)(p::Perm{I}) where {I} =
    Perm(I[ihom[f^p] for f in ihom.features])

function (ihom::InducingHomomorphism)(orb::O) where {T,O<:AbstractOrbit{T,Nothing}}
    elts = ihom.(orb)
    dict = Dict(e => nothing for e in elts)
    return O(elts, dict)
end

function (ihom::InducingHomomorphism)(χ::CF) where {CF<:ClassFunction}
    iccG = ihom.(conjugacy_classes(χ))
    return CF(values(χ), χ.inv_of, iccG)
end
