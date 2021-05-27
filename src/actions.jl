abstract type Action end
abstract type ByPermutations <: Action end
abstract type ByLinearTransformation <: Action end

abstract type InducedActionHomomorphism{T} end

function (hom::InducedActionHomomorphism)(orb::O) where {T,O<:AbstractOrbit{T,Nothing}}
    elts = hom.(orb)
    dict = Dict(e => nothing for e in elts)
    return Orbit(elts, dict)
end

function (hom::InducedActionHomomorphism)(chars::AbstractArray{<:Character})
    χ = first(chars)
    iccG = hom.(conjugacy_classes(χ))
    return [Character(values(χ), χ.inv_of, iccG) for χ in chars]
end

function induce(ac::Action, hom::InducedActionHomomorphism, g::GroupElement)
    throw(
        "No fallback is provided for $(typeof(ac)). You need to implement `induce(::$(typeof(ac)), ::$(typeof(hom)), ::$(typeof(g)))`.",
    )
end

decompose(x::Any, hom::InducedActionHomomorphism) = throw(
    "No fallback is provided for $(typeof(x)). You need to implement `decompose(::$(typeof(x)), ::$(typeof(hom)))`.",
)

coeff_type(ac::ByLinearTransformation) = throw(
    "No fallback is provided for $(typeof(ac)). You need to implement `coeff_type(::$(typeof(ac)))`.",
)

#=
implements:
* features(action_hom)
* action(hom::InducedActionHomomorphism) → Action
* action(hom::InducedActionHomomorphism, g::GroupElement, feature) → return action of `g` on `feature`
=#

struct ExtensionHomomorphism{T,A,I,V} <: InducedActionHomomorphism{T}
    ac::A
    features::V # supports linear indexing
    reversef::Dict{T,I}
end

if VERSION < v"1.3.0"
    hom(::ExtensionHomomorphism)(g::GroupElement) = induce(action(hom), hom, g)
else
    (hom::InducedActionHomomorphism)(g::GroupElement) = induce(action(hom), hom, g)
end

function ExtensionHomomorphism{I}(ac::Action, features) where {I<:Integer}
    @assert typemax(I) >= length(features) "Use wider Integer type!"
    reversef = Dict{eltype(features),I}(f => idx for (idx, f) in pairs(features))
    return ExtensionHomomorphism(ac, features, reversef)
end

ExtensionHomomorphism(ac::Action, features) = ExtensionHomomorphism{Int16}(ac, features)

# interface:
features(hom::ExtensionHomomorphism) = hom.features
action(hom::ExtensionHomomorphism) = hom.ac

# user convenience functions:
Base.getindex(ehom::ExtensionHomomorphism, i::Integer) = ehom.features[i]
Base.getindex(ehom::ExtensionHomomorphism{T}, f::T) where {T} = ehom.reversef[f]

function induce(::ByPermutations, hom::ExtensionHomomorphism, g::GroupElement)
    return Perm(vec([hom[action(action(hom), g, f)] for f in features(hom)]))
end

using SparseArrays
function induce(ac::ByLinearTransformation, hom::InducedActionHomomorphism, g::GroupElement)

    n = length(features(hom))

    I = Int[]
    J = Int[]
    V = coeff_type(ac)[]

    for (i, f) in enumerate(features(hom))
        k = action(action(hom), g, f)
        idcs, vals = decompose(k, hom)
        append!(I, fill(i, length(idcs)))
        append!(J, idcs)
        append!(V, vals)
    end
    return sparse(I, J, V)
end

# ala permutation action
function decompose(m::T, hom::ExtensionHomomorphism{T}) where {T}
    return [hom[m]], [one(coeff_type(hom.ac))]
end
