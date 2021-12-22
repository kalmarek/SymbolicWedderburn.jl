abstract type Action end
abstract type ByPermutations <: Action end
abstract type ByLinearTransformation <: Action end

abstract type InducedActionHomomorphism{A,T} end
#=
implements:
* basis(action_hom)
* action(hom::InducedActionHomomorphism) → Action
* action(hom::InducedActionHomomorphism, g::GroupElement, feature) → the action of `g` on `feature`
=#

function action(
    hom::InducedActionHomomorphism{<:ByPermutations},
    g::GroupElement,
    v::AbstractVector,
)
    return v^induce(hom, g)
end

function action(
    hom::InducedActionHomomorphism{<:ByLinearTransformation},
    g::GroupElement,
    v::AbstractVector,
)
    return induce(hom, g) * v
end

function induce(ac::Action, hom::InducedActionHomomorphism, g::GroupElement)
    return throw(
        "No fallback is provided for $(typeof(ac)). You need to implement `induce(::$(typeof(ac)), ::$(typeof(hom)), ::$(typeof(g)))`.",
    )
end

induce(hom::InducedActionHomomorphism, g::GroupElement) =
    induce(action(hom), hom, g)

decompose(x::Any, hom::InducedActionHomomorphism) = throw(
    "No fallback is provided for $(typeof(x)). You need to implement `decompose(::$(typeof(x)), ::$(typeof(hom)))`.",
)
coeff_type(ac::ByLinearTransformation) = throw(
    "No fallback is provided for $(typeof(ac)). You need to implement `_coeff_type(::$(typeof(ac)))`.",
)

# Exceeding typemax(UInt16) here would mean e.g. that you're trying to block-diagonalize
# an SDP constraint of size 65535×65535, which is highly unlikely ;)
_int_type(::Type{<:InducedActionHomomorphism}) = UInt16
_int_type(hom::InducedActionHomomorphism) = _int_type(typeof(hom))

function induce(
    ac::ByLinearTransformation,
    hom::InducedActionHomomorphism,
    g::GroupElement,
)
    n = length(basis(hom))

    I = Int[]
    J = Int[]
    V = coeff_type(ac)[]

    for (i, f) in enumerate(basis(hom))
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
    return [hom[m]], [one(coeff_type(action(hom)))]
end

struct ExtensionHomomorphism{A<:Action,T,B<:StarAlgebras.AbstractBasis{T}} <:
       InducedActionHomomorphism{A,T}
    action::A
    basis::B
end

# see the comment about UInt16 above
ExtensionHomomorphism(action::Action, basis) =
    ExtensionHomomorphism(_int_type(ExtensionHomomorphism), action, basis)

ExtensionHomomorphism(::Type{I}, action::Action, basis) where {I<:Integer} =
    ExtensionHomomorphism(action, StarAlgebras.Basis{I}(basis))

_int_type(::Type{<:ExtensionHomomorphism{A,T,B}}) where {A,T,B} = _int_type(B)
_int_type(::Type{<:StarAlgebras.AbstractBasis{T,I}}) where {T,I} = I

Base.getindex(ehom::ExtensionHomomorphism, i::Integer) = basis(ehom)[i]
Base.getindex(ehom::ExtensionHomomorphism{A,T}, f::T) where {A,T} =
    ehom.basis[f]

# interface:
StarAlgebras.basis(hom::ExtensionHomomorphism) = hom.basis
action(hom::ExtensionHomomorphism) = hom.action

function induce(::ByPermutations, hom::ExtensionHomomorphism, g::GroupElement)
    I = _int_type(hom)
    return Perm(vec(I[hom[action(action(hom), g, f)] for f in basis(hom)]))
end

struct CachedExtensionHomomorphism{A,T,G,H,E<:InducedActionHomomorphism{A,T}} <:
       InducedActionHomomorphism{A,T}
    ehom::E
    cache::Dict{G,H}
end

action(h::CachedExtensionHomomorphism) = action(h.ehom)
StarAlgebras.basis(h::CachedExtensionHomomorphism) = basis(h.ehom)
_int_type(chom::CachedExtensionHomomorphism) = _int_type(h.ehom)

Base.getindex(h::CachedExtensionHomomorphism, x::Any) = h.ehom[x]

function CachedExtensionHomomorphism(
    G::Group,
    action::Action,
    basis;
    precompute = false,
)
    hom = ExtensionHomomorphism(action, basis)
    S = typeof(induce(hom, one(G)))
    chom = CachedExtensionHomomorphism{eltype(G),S}(hom)
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

function _induce(ac::Action, hom::CachedExtensionHomomorphism, g::GroupElement)
    if !haskey(hom.cache, g)
        hom.cache[g] = induce(ac, hom.ehom, g)
    end
    return hom.cache[g]
end

induce(ac::Action, hom::CachedExtensionHomomorphism, g::GroupElement) =
    _induce(ac, hom, g)
# disabmiguation:
induce(
    ac::ByLinearTransformation,
    hom::CachedExtensionHomomorphism,
    g::GroupElement,
) = _induce(ac, hom, g)
