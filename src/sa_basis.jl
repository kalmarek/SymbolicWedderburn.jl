"""
    isotypical_basis(χ::AbstractClassFunction)
Compute a basis of the image of the projection corresponding to a class function `χ`.

Return the coefficients of basis vectors in an invariant subspace corresponding to `χ`
(so called _isotypical subspace_) in the action on `ℝⁿ` encoded by the conjugacy classes of `χ`.
"""
function isotypical_basis(χ::AbstractClassFunction)
    u, weight = matrix_projection(χ)
    image, pivots = if iszero(weight)
        u_ = similar(u, 0, size(u, 2))
        image_basis!(u_)
    else
        u .*= weight
        image_basis!(u)
    end
    @debug "isotypical subspace corresponding to χ has dimension $(length(pivots))" χ

    return image[pivots, :]
end

struct SemisimpleSummand{T}
    basis::T
    multiplicity::Int
end

basis(b::SemisimpleSummand) = b.basis
Base.convert(::Type{M}, b::SemisimpleSummand) where {M<:AbstractMatrix} =
    convert(M, basis(b))

multiplicity(b::SemisimpleSummand) = b.multiplicity

function degree(b::SemisimpleSummand)
    d, r = divrem(size(basis(b), 1), multiplicity(b))
    @assert iszero(r)
    return d
end

function extended_characters(::Type{T}, G::Group, basis, action) where {T}
    chars = characters_dixon(T, G)
    ehom = ExtensionHomomorphism(action, basis)
    return ehom(chars) # result: Cyclotomic{T} valued Characters
end

"""
    symmetry_adapted_basis([T::Type,] G::AbstractPermutationGroup[, S=Rational{Int}])
Compute a basis for the linear space `ℝⁿ` which is invariant under the symmetry of `G`.

The permutation group is acting naturally on `1:degree(G)`. The coefficients of
the invariant basis are returned in (orthogonal) blocks corresponding to irreducible
characters of `G`.

Arguments:
* `S` controlls the types of `Cyclotomic`s used in the computation of
character table. Exact type are preferred. For larger groups `G` `Rational{BigInt}`
might be necessary.
* `T` controls the type of coefficients of the returned basis.

!!! Note:
Each block is invariant under the action of `G`, i.e. the action may permute
vectors from symmetry adapted basis within each block. The blocks are guaranteed
to be orthogonal. If `T<:LinearAlgebra.BlasFloat` BLAS routines will be used to
orthogonalize vectors within each block.
"""
symmetry_adapted_basis(G::AbstractPermutationGroup, S::Type = Rational{Int}) =
    symmetry_adapted_basis(characters_dixon(S, G))

symmetry_adapted_basis(T::Type, G::AbstractPermutationGroup, S::Type = Rational{Int}) =
    symmetry_adapted_basis(T, characters_dixon(S, G))

"""
    symmetry_adapted_basis([T::Type,] G::Group, basis, action[, S=Rational{Int}])
Compute a basis for the linear space spanned by `basis` which is invariant under
the symmetry of `G`.

* The action used in these computations is
> `(b,g) → action(b,g)` for `b ∈ basis`, `g ∈ G`
and needs to be defined by the user.
* It is assumed that `G` acts on a subset of basis and the action needs to be
extended to the whole `basis`. If `G` is a permutation group already acting on
the whole `basis`, a call to `symmetry_adapted_basis(G)` is preferred.
* For inducing the action `basis` needs to be indexable and iterable
(e.g. in the form of an `AbstractVector`).

Arguments:
* `S` controlls the types of `Cyclotomic`s used in the computation of
character table. Exact type are preferred. For larger groups `G` `Rational{BigInt}`
might be necessary.
* `T` controls the type of coefficients of the returned basis.

!!! Note:
Each block is invariant under the action of `G`, i.e. the action may permute
vectors from symmetry adapted basis within each block. The blocks are guaranteed
to be orthogonal. If `T<:LinearAlgebra.BlasFloat` BLAS routines will be used to
orthogonalize vectors within each block.
"""
function symmetry_adapted_basis(G::Group, basis, action, S::Type = Rational{Int})
    chars_ext = extended_characters(S, G, basis, action)
    return symmetry_adapted_basis(chars_ext)
end

function symmetry_adapted_basis(
    ::Type{T},
    G::Group,
    basis,
    action,
    ::Type{S} = Rational{Int},
) where {T,S}
    chars_ext = extended_characters(S, G, basis, action)
    return symmetry_adapted_basis(T, chars_ext)
end

function symmetry_adapted_basis(
    ::Type{T},
    chars::AbstractVector{<:AbstractClassFunction},
) where {T}
    chars = T <: Real ? affordable_real!(deepcopy.(chars)) : chars
    return symmetry_adapted_basis(Character{T}.(chars))
end

_multiplicities(ψ, chars) = Int[Int(dot(ψ, χ)) / Int(dot(χ, χ)) for χ in chars]
_multiplicities(ψ, chars::AbstractVector{<:AbstractClassFunction{T}}) where T<:AbstractFloat =
    [round(Int, dot(ψ, χ) / dot(χ, χ)) for χ in chars]
_multiplicities(ψ, chars::AbstractVector{<:AbstractClassFunction{T}}) where T<:Complex =
    [round(Int, real(dot(ψ, χ) / dot(χ, χ))) for χ in chars]

function symmetry_adapted_basis(chars::AbstractVector{<:AbstractClassFunction})
    ψ = let χ = first(chars)
        action_character(conjugacy_classes(χ), χ.inv_of)
    end

    multiplicities = _multiplicities(ψ, chars)
    degrees = degree.(chars)

    @debug "Decomposition into character spaces:
    degrees:        $(join([lpad(d, 6) for d in degrees], ""))
    multiplicities: $(join([lpad(m, 6) for m in multiplicities], ""))"

    dot(multiplicities, degrees) == degree(ψ) ||
        @error "characters do not constitute a complete basis for action:
        $(dot(multiplicities, degrees)) ≠ $(degree(ψ))"

    args = filter(!iszero ∘ last, collect(zip(chars, multiplicities)))

    res = map(args) do (χ, m)
        Threads.@spawn SemisimpleSummand(isotypical_basis(χ), m)
    end

    return fetch.(res)
end
