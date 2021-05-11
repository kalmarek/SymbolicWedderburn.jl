abstract type AbstractActionHomomorphism{T} end

struct ExtensionHomomorphism{T,I<:Integer,V,Op} <: AbstractActionHomomorphism{T}
    features::V
    reversef::Dict{T,I}
    op::Op
end

function ExtensionHomomorphism{I}(features, op) where I<:Integer
    @assert typemax(I) >= length(features) "Use wider Integer type!"
    reversef = Dict{eltype(features), I}(f => idx for (idx, f) in pairs(features))
    return ExtensionHomomorphism(features, reversef, op)
end

ExtensionHomomorphism(features, op) = ExtensionHomomorphism{Int16}(features, op)

Base.getindex(ehom::ExtensionHomomorphism, i::Integer) = ehom.features[i]
Base.getindex(ehom::ExtensionHomomorphism{T}, f::T) where {T} = ehom.reversef[f]

function permutation(ehom::SymbolicWedderburn.ExtensionHomomorphism{T}, els::Vector{T}) where {T}
    return Perm([ehom[el] for el in els])
end

function (ehom::SymbolicWedderburn.ExtensionHomomorphism)(g::GroupsCore.GroupElement)
    return permutation(ehom, [ehom.op(f, g) for f in ehom.features])
end

function (ehom::ExtensionHomomorphism)(orb::O) where {T,O<:AbstractOrbit{T,Nothing}}
    elts = ehom.(orb)
    dict = Dict(e => nothing for e in elts)
    return Orbit(elts, dict)
end

function (ehom::ExtensionHomomorphism)(χ::Character)
    iccG = ehom.(conjugacy_classes(χ))
    return Character(values(χ), χ.inv_of, iccG)
end
