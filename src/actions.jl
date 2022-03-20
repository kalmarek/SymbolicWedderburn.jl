abstract type Action end

"""
    ByPermutations <: Action
A type of action where a group acts through permutations on a set.

It means that for every group element `g ∈ G` and every `s ∈ S`
```
action(act::ByPermutations, g, s) = s′
```
where `s′ ∈ S` is a (potentially different) element of `S`. By fixing an order
of elements in `S` one can then represent `g` as a permutation of degree `|S|`.
"""
abstract type ByPermutations <: Action end

"""
    ByLinearTransformation <: Action
A type of action where a group acts through linear transformations on an
(implicit) Euclidean space `ℝᴮ` with basis `B = (e₁, … eₙ)`.

That means that for every group element `g ∈ G` and every `eₖ ∈ B`
```
action(act::ByPermutations, g, eₖ) = v
```
where `v = a₁e₁ + ⋯ + aₙeₙ`. Additionally [`decompose`](@ref) must be
implemented to return a (possibly sparse) decomposition of `v` in `B`.
"""
abstract type ByLinearTransformation <: Action end

"""
    BySignedPermutations
A type of action where a group acts through permutations _with sign_ on an
(implicit) Euclidean space `ℝᴮ` with basis `B = (e₁, …, eₙ)`.

It means that for every group element `g ∈ G` and every `eₖ ∈ B` the result of
action of `g` on `eₖ` is `u·eₗ` for some `1≤l≤n`, where `u` is a root of unity
(i.e. `±1` in the real case). To accomplish this it is required that
```
action(act::BySignedPermutations, g, eₖ) == (eₗ, u)
```

!!! warning
    Only support for the real case (i.e. `u = ±1`) is implemented at the moment.

"""
abstract type BySignedPermutations <: ByLinearTransformation end

abstract type InducedActionHomomorphism{A,T} end

#=
implements:
* basis(action_hom)
* action(hom::InducedActionHomomorphism) → Action
* action(hom::InducedActionHomomorphism, g::GroupElement, feature) → the action of `g` on `feature`
=#

"""
    action(hom::InducedActionHomomorphism, g::GroupElement, x)
Return the result of `g` acting on `x` through action homomorphism `hom`.

This can be understood as first evaluating the homomorphism: `h = hom(g)` and
then computing `x^h`, the action of the result on `x`.
"""
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
        "No fallback is provided for $(typeof(ac)).
        You need to implement
        `induce(::$(typeof(ac)), ::$(typeof(hom)), ::$(typeof(g)))`.",
    )
end

induce(hom::InducedActionHomomorphism, g::GroupElement) =
    induce(action(hom), hom, g)

"""
    decompose(v, hom::InducedActionHomomorphism)
Decompose element `v` in the (implicit) basis provided by `hom`.

Let `B = basis(hom)::StarAlgebras.AbstractBasis`. Then `v` should be decomposed
as a unique linear combination `v = a₁b₁ + ⋯ aₙbₙ` and the indices of `bᵢ`s in `B`
and a vector of coefficients `A` returned, so that

```
@assert sparsevec(decompose(v, hom)...) == v
```

!!! note
    For performance reasons it is best to drop zeros in `A`, i.e. return a
    sparse representation.

See also [`ByLinearTransformation`](@ref).
"""
decompose(x::Any, hom::InducedActionHomomorphism) = throw(
    "No fallback is provided for $(typeof(x)). You need to implement
    `decompose(::$(typeof(x)), ::$(typeof(hom)))`.",
)
coeff_type(ac::ByLinearTransformation) = throw(
    "No fallback is provided for $(typeof(ac)). You need to implement
    `coeff_type(::$(typeof(ac)))`.",
)
coeff_type(hom::InducedActionHomomorphism) = coeff_type(action(hom))
coeff_type(::ByPermutations) = Int

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
    lock::Base.Threads.SpinLock
end

CachedExtensionHomomorphism{G,H}(hom::InducedActionHomomorphism) where {G,H} =
    CachedExtensionHomomorphism(hom, Dict{G,H}(), Threads.SpinLock())

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
    @sync if precompute
        for g in G
            Threads.@spawn begin
                induce(action, chom, g)
            end
        end
    end
    return chom
end

function _induce(ac::Action, chom::CachedExtensionHomomorphism, g::GroupElement)
    if !haskey(chom.cache, g)
        val = induce(ac, chom.ehom, g)
        lock(chom.lock) do
            chom.cache[g] = val
        end
    end
    return chom.cache[g]
end

induce(ac::Action, hom::CachedExtensionHomomorphism, g::GroupElement) =
    _induce(ac, hom, g)
# disabmiguation:
induce(
    ac::ByLinearTransformation,
    hom::CachedExtensionHomomorphism,
    g::GroupElement,
) = _induce(ac, hom, g)

function _induce(
    ac::BySignedPermutations,
    hom::InducedActionHomomorphism,
    g::GroupElement, # need to be signed permutation group element
)
    bs = basis(hom)
    n = length(bs)

    I = Int[]
    J = Int[]
    V = coeff_type(ac)[]

    for (i, f) in enumerate(bs)
        k, s = action(action(hom), g, f) # k::MonoidElement
        push!(I, i)
        push!(J, bs[k]) # refactor into decompose(::MonoidElement, hom)?
        push!(V, s)
    end
    return sparse(I, J, V)
end

# disambiguation:
induce(
    ac::BySignedPermutations,
    hom::CachedExtensionHomomorphism,
    g::GroupElement,
) = _induce(ac, hom, g)

coeff_type(::BySignedPermutations) = Int # lets not worry about roots of unity
