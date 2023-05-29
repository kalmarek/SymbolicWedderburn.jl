# abstract type Action is defined in ext_homomorphisms

"""
    GroupActionError
`GroupActionError` is an indicator that the group action defined by
`(G, action, basis)` might be incorrect.

If you implement the action in question yourself or encounter inconsistent results
is advisable to run `check_group_action(G, action, basis, full_check=true)` to
check the correctness of the implementation of the action.
"""
struct GroupActionError <: Exception
    msg::String
end

function Base.showerror(io::IO, err::GroupActionError)
    print(io, "Failure of group action axiom:")
    return print(io, err.msg)
end

function __check_group_action_axioms(itr, act::Action, set_with_action)
    id = one(first(itr))
    for x in set_with_action
        if action(act, id, x) != x
            throw(
                GroupActionError("`one(G)` doesn't act as identity on `$(x)`"),
            )
        end
        for g in itr, h in itr
            a = action(act, g * h, x)
            b = action(act, h, action(act, g, x))
            if a != b
                throw(
                    GroupActionError(
                        "action by `g=$(g)` and `h=$(h)` fail right-associativity on `x=$(x)` :\n " *
                        "`action(act, g * h, x) = $a ≠ $b = action(act, h, action(act, g, x))`",
                    ),
                )
            end
        end
    end
    return true
end

function check_group_action(
    G::Group,
    act::Action,
    basis::AbstractVector;
    full_check = false,
    elts = full_check ? (1:length(basis)) : (1:min(length(basis), ngens(G))),
)
    X = @view basis[elts]
    if full_check
        return __check_group_action_axioms(G, act, X)
    else
        return __check_group_action_axioms(gens(G), act, X)
    end
end

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
(implicit) Euclidean space `ℝᴮ` with basis `B = (e₁, … eₙ)`.

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

function coeff_type(ac::ByLinearTransformation)
    throw("No fallback is provided for $(typeof(ac)). You need to implement
          `coeff_type(::$(typeof(ac)))`.")
end

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
function induce(
    ac::ByLinearTransformation,
    hom::CachedExtensionHomomorphism,
    g::GroupElement,
)
    return _induce(ac, hom, g)
end

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
function decompose(x, hom::InducedActionHomomorphism)
    throw("""No fallback is provided for $(typeof(x)). You need to implement
          `decompose(::$(typeof(x)), ::$(typeof(hom)))`.""")
end

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

    Please consult the docstring of [ByLinearTransformation](@ref) for the necessary methods.
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
        k, s = action(action(hom), g, f)
        push!(I, i)
        push!(J, bs[k])
        push!(V, s)
    end
    return sparse(I, J, V)
end

# disabmiguation
function induce(
    ac::BySignedPermutations,
    hom::CachedExtensionHomomorphism,
    g::GroupElement,
)
    return _induce(ac, hom, g)
end
