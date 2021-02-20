abstract type AbstractActionHomomorphism{T} end

struct ExtensionHomomorphism{T,V,Op} <: AbstractActionHomomorphism{T}
    features::V
    reversef::Dict{T,Int}
    op::Op
end

ExtensionHomomorphism(features, op) = ExtensionHomomorphism(
    features,
    Dict(f => idx for (idx, f) in enumerate(features)),
    op,
)

Base.getindex(ehom::ExtensionHomomorphism, i::Integer) = ehom.features[i]
Base.getindex(ehom::ExtensionHomomorphism{T}, f::T) where {T} = ehom.reversef[f]
(ehom::ExtensionHomomorphism)(p::Perm{I}) where {I} =
    Perm(vec(I[ehom[ehom.op(f, p)] for f in ehom.features]))

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
