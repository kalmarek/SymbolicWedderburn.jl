# abstract type Action is defined in ext_homomorphisms

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

coeff_type(::ByPermutations) = Int

function induce(::ByPermutations, hom::ExtensionHomomorphism, g::GroupElement)
    I = _int_type(hom)
    return Perm(vec(I[hom[action(action(hom), g, f)] for f in basis(hom)]))
end

"""
    ByLinearTransformation <: Action
A type of action where a group acts through linear transformations on an
(implicit) linear space `ℝᴮ` with basis `B = (e₁, … eₙ)`.

That means that for every group element `g ∈ G` and every `eₖ ∈ B`
```
action(act::ByLinearTransformation, g, eₖ) = v
```
where `v = a₁e₁ + ⋯ + aₙeₙ`. Additionally [`decompose`](@ref) must be
implemented to return a (possibly sparse) decomposition of `v` in `B`.
"""
abstract type ByLinearTransformation <: Action end

function action(
    hom::InducedActionHomomorphism{<:ByLinearTransformation},
    g::GroupElement,
    v::AbstractVector,
)
    return induce(hom, g) * v
end

coeff_type(ac::ByLinearTransformation) =
    throw("No fallback is provided for $(typeof(ac)). You need to implement
          `coeff_type(::$(typeof(ac)))`.")

"""
    induce(ac::ByLinearTransformation, hom::InducedActionHomomorphism, g::GroupElement)
Return a sparse matrix representation of the action by `g`` in the basis provided by `hom`
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

# disabmiguation:
induce(
    ac::ByLinearTransformation,
    hom::CachedExtensionHomomorphism,
    g::GroupElement,
) = _induce(ac, hom, g)

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
decompose(x, hom::InducedActionHomomorphism) =
    throw("""No fallback is provided for $(typeof(x)). You need to implement
          `decompose(::$(typeof(x)), ::$(typeof(hom)))`.""")

"""
    BySignedPermutations <: ByLinearTransformation
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

coeff_type(::BySignedPermutations) = Int # lets not worry about roots of unity

function induce(
    ac::BySignedPermutations,
    hom::InducedActionHomomorphism,
    g::GroupElement,
)
    bs = basis(hom)
    n = length(bs)

    I = Int[]
    J = Int[]
    V = coeff_type(ac)[]

    for (i, f) in enumerate(bs)
        k, s = action(action(hom), g, f) # act and decompose
        push!(I, i)
        push!(J, bs[k])
        push!(V, s)
    end
    return sparse(I, J, V)
end

# disabmiguation
induce(
    ac::BySignedPermutations,
    hom::CachedExtensionHomomorphism,
    g::GroupElement,
) = _induce(ac, hom, g)