abstract type Action end
abstract type ByPermutations <: Action end
abstract type ByLinearTransformation <: Action end

abstract type InducedActionHomomorphism{A,T} end
#=
implements:
* features(action_hom)
* action(hom::InducedActionHomomorphism) → Action
* action(hom::InducedActionHomomorphism, g::GroupElement, feature) → the action of `g` on `feature`
=#

function induce(ac::Action, hom::InducedActionHomomorphism, g::GroupElement)
    throw(
        "No fallback is provided for $(typeof(ac)). You need to implement `induce(::$(typeof(ac)), ::$(typeof(hom)), ::$(typeof(g)))`.",
    )
end

induce(hom::InducedActionHomomorphism, g::GroupElement) = induce(action(hom), hom, g)

decompose(x::Any, hom::InducedActionHomomorphism) = throw(
    "No fallback is provided for $(typeof(x)). You need to implement `decompose(::$(typeof(x)), ::$(typeof(hom)))`.",
)

_coeff_type(ac::ByLinearTransformation) = throw(
    "No fallback is provided for $(typeof(ac)). You need to implement `_coeff_type(::$(typeof(ac)))`.",
)

# Exceeding typemax(UInt16) here would mean e.g. that you're trying to block-diagonalize
# an SDP constraint of size 65535×65535, which is highly unlikely ;)
_int_type(::InducedActionHomomorphism) = UInt16

function induce(ac::ByLinearTransformation, hom::InducedActionHomomorphism, g::GroupElement)

    n = length(features(hom))

    I = Int[]
    J = Int[]
    V = _coeff_type(ac)[]

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
function decompose(m::T, hom::InducedActionHomomorphism{A,T}) where {A,T}
    return [hom[m]], [one(_coeff_type(action(hom)))]
end


struct ExtensionHomomorphism{A,T,I,V} <: InducedActionHomomorphism{A,T}
    ac::A
    features::V # supports linear indexing
    reversef::Dict{T,I}
end

# see the comment about UInt16 above
ExtensionHomomorphism(ac::Action, features) = ExtensionHomomorphism{UInt16}(ac, features)

function ExtensionHomomorphism{I}(ac::Action, features) where {I <: Integer}
    @assert typemax(I) >= length(features) "Use wider Integer type!"
    reversef = Dict{eltype(features),I}(f => idx for (idx, f) in pairs(features))
    return ExtensionHomomorphism(ac, features, reversef)
end

_int_type(::ExtensionHomomorphism{A,T,I}) where {A,T,I} = I

Base.getindex(ehom::ExtensionHomomorphism, i::Integer) = ehom.features[i]
Base.getindex(ehom::ExtensionHomomorphism{A,T}, f::T) where {A,T} = ehom.reversef[f]

# interface:
features(hom::ExtensionHomomorphism) = hom.features
action(hom::ExtensionHomomorphism) = hom.ac

function induce(::ByPermutations, hom::ExtensionHomomorphism, g::GroupElement)
    I = _int_type(hom)
    return Perm(vec(I[hom[action(action(hom), g, f)] for f in features(hom)]))
end

struct CachedExtensionHomomorphism{A,T,G,H,E <: InducedActionHomomorphism{A,T}} <: InducedActionHomomorphism{A,T}
    ehom::E
    cache::Dict{G,H}
end

for f in (:_int_type, :features, :action)
    @eval $f(h::CachedExtensionHomomorphism) = $f(h.ehom)
end

function CachedExtensionHomomorphism(G::Group, action::Action, basis; precompute=false)
    hom = ExtensionHomomorphism(action, basis)
    S = typeof(induce(hom, one(G)))
    chom = CachedExtensionHomomorphism{eltype(G), S}(hom)
    if precompute
        # to make sure that the future access to chom is read-only, i.e. thread-safe
        # one may choose to precompute the images
        for g in G
            induce(chom, g)
        end
    end
    return chom
end

CachedExtensionHomomorphism{G,H}(hom::InducedActionHomomorphism) where {G,H} =
    CachedExtensionHomomorphism(hom, Dict{G,H}())

function induce(ac::Action, hom::CachedExtensionHomomorphism, g::GroupElement)
    if !haskey(hom.cache, g)
        hom.cache[g] = induce(ac, hom.ehom, g)
    end
    return hom.cache[g]
end
