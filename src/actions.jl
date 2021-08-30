abstract type Action end
abstract type ByPermutations <: Action end
abstract type ByLinearTransformation <: Action end

abstract type InducedActionHomomorphism{A, T} end

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

struct ExtensionHomomorphism{A,T,I,V} <: InducedActionHomomorphism{A, T}
    ac::A
    features::V # supports linear indexing
    reversef::Dict{T,I}
end

_inttype(::ExtensionHomomorphism{A,T,I}) where {A,T,I} = I

@static if VERSION < v"1.3.0"
    # cannot add methods to an abstract type in pre Julia v1.3
    (hom::ExtensionHomomorphism)(g::GroupElement) = induce(action(hom), hom, g)

    function (hom::ExtensionHomomorphism)(orb::O) where {T,O<:AbstractOrbit{T,Nothing}}
        elts = hom.(orb)
        dict = Dict(e => nothing for e in elts)
        return Orbit(elts, dict)
    end

    function (hom::ExtensionHomomorphism)(chars::AbstractArray{<:Character})
        χ = first(chars)
        iccG = hom.(conjugacy_classes(χ))
        return [Character(values(χ), χ.inv_of, iccG) for χ in chars]
    end
else
    (hom::InducedActionHomomorphism)(g::GroupElement) = induce(action(hom), hom, g)

    function (hom::InducedActionHomomorphism)(orb::O) where {T,O<:AbstractOrbit{T,Nothing}}
        S = typeof(hom(one(first(orb))))
        elts = Vector{S}(undef, length(orb))
        for i in 1:length(orb)
            elts[i] = hom(orb.elts[i])
        end
        return Orbit(elts)
    end

    function (hom::InducedActionHomomorphism)(chars::AbstractArray{<:Character})
        χ = first(chars)
        iccG = hom.(conjugacy_classes(χ))
        return [Character(values(χ), χ.inv_of, iccG) for χ in chars]
    end
end

function ExtensionHomomorphism{I}(ac::Action, features) where {I<:Integer}
    @assert typemax(I) >= length(features) "Use wider Integer type!"
    reversef = Dict{eltype(features),I}(f => idx for (idx, f) in pairs(features))
    return ExtensionHomomorphism(ac, features, reversef)
end

# Exceeding typemax(UInt16) here would mean e.g. that you're trying to block-diagonalize SDP constraint of size 65535×65535, which is highly unlikely ;)
ExtensionHomomorphism(ac::Action, features) = ExtensionHomomorphism{UInt16}(ac, features)

# interface:
features(hom::ExtensionHomomorphism) = hom.features
action(hom::ExtensionHomomorphism) = hom.ac

# user convenience functions:
Base.getindex(ehom::ExtensionHomomorphism, i::Integer) = ehom.features[i]
Base.getindex(ehom::ExtensionHomomorphism{A,T}, f::T) where {A,T} = ehom.reversef[f]

function induce(::ByPermutations, hom::ExtensionHomomorphism, g::GroupElement)
    I = _inttype(hom)
    return Perm(vec(I[hom[action(action(hom), g, f)] for f in features(hom)]))
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
function decompose(m::T, hom::InducedActionHomomorphism{A, T}) where {A, T}
    return [hom[m]], [one(coeff_type(action(hom)))]
end

struct CachedExtensionHomomorphism{A, T, G, H, E<:InducedActionHomomorphism{A, T}} <: InducedActionHomomorphism{A,T}
    ehom::E
    cache::Dict{G, H}
end

for f in (:_inttype, :features, :action)
    @eval $f(h::CachedExtensionHomomorphism) = $f(h.ehom)
end

CachedExtensionHomomorphism{G, H}(hom::InducedActionHomomorphism) where {G,H} =
    CachedExtensionHomomorphism(hom, Dict{G, H}())

function induce(ac::Action, hom::CachedExtensionHomomorphism, g::GroupElement)
    if !haskey(hom.cache, g)
        hom.cache[g] = induce(ac, hom.ehom, g)
    end
    return hom.cache[g]
end
