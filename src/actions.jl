abstract type AbstractActionHomomorphism{T} end

struct ExtensionHomomorphism{T,V} <: AbstractActionHomomorphism{T}
    features::V
    reversef::Dict{T,Int}
end

ExtensionHomomorphism(features) = ExtensionHomomorphism(
    features,
    Dict(f => idx for (idx, f) in enumerate(features)),
)

Base.getindex(ehom::ExtensionHomomorphism, i::Integer) = ehom.features[i]
Base.getindex(ehom::ExtensionHomomorphism{T}, f::T) where {T} = ehom.reversef[f]
(ehom::ExtensionHomomorphism)(p::Perm{I}) where {I} =
    Perm(vec(I[ehom[f^p] for f in ehom.features]))

function (ehom::ExtensionHomomorphism)(
    orb::O,
) where {T,O<:AbstractOrbit{T,Nothing}}
    elts = ehom.(orb)
    dict = Dict(e => nothing for e in elts)
    return O(elts, dict)
end

function (ehom::ExtensionHomomorphism)(χ::CF) where {CF<:ClassFunction}
    iccG = ehom.(conjugacy_classes(χ))
    return CF(values(χ), χ.inv_of, iccG)
end
