abstract type Action end
abstract type ByPermutations <: Action end
abstract type ByLinearTransformation <: Action end
abstract type BySignedPermutations <: ByLinearTransformation end
abstract type BySigns <: ByLinearTransformation end
abstract type InducedActionHomomorphism{A,T} end

"""
    BySignedPermutations

Instances of `BySignedPermutations` must satisfy
```
action(act::BySignedPermutations, g, b) == (gb, s)
```
whenever the action maps `(g, b)` to  `s*gb`. Note: `s` must be a sign (this limitation is a design choice).

For problems over the complex field, `s` is necessarily a root of unity. You need to define your own abstract type.
"""

#=
implements:
* basis(action_hom)
* action(hom::InducedActionHomomorphism) → Action
* action(hom::InducedActionHomomorphism, g::GroupElement, feature) → the action of `g` on `feature`
=#

"""
Compute the vector g⋅v which is the result of acting g on v by permutation.
"""
function action(
    hom::InducedActionHomomorphism{<:ByPermutations},
    g::GroupElement,
    v::AbstractVector,
)
    return v^induce(hom, g)
end

"""
Compute the vector g⋅v which is the result of acting g on v by linear transformation.
"""
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
coeff_type(hom::InducedActionHomomorphism) = coeff_type(action(hom))
coeff_type(::ByPermutations) = Int

# Exceeding typemax(UInt16) here would mean e.g. that you're trying to block-diagonalize
# an SDP constraint of size 65535×65535, which is highly unlikely ;)
_int_type(::Type{<:InducedActionHomomorphism}) = UInt16
_int_type(hom::InducedActionHomomorphism) = _int_type(typeof(hom))

"""
Compute the matrix corresponding to g acting by linear transformation.
Note: the return matrix is of type SparseMatrixCSC but may in fact be dense.
"""
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
        lck = Threads.SpinLock()
        tasks = [
        Threads.@spawn begin
            val = SymbolicWedderburn.induce(action, chom.ehom, g)
            lock(lck) do
                chom.cache[g] = val
            end
        end for g in G]
        foreach(wait, tasks)
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


"""
Compute the vector g⋅v which is the result of acting g on v by signed permutation.
"""
function action(
    hom::InducedActionHomomorphism{<:BySignedPermutations},
    g::GroupElement,# TODO: narrow to signed permutation group element
    v::AbstractVector,
)
    return action(ByPermutations(), g, v), action(BySigns(), g, v)
end

# TODO: define action BySigns, which requires API for signed permutation group


"""
Compute the sparse matrix representing the action of g by signed permutation.
"""
function induce(
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

coeff_type(::BySigns) = Int
