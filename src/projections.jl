"""
    matrix_projection(χ::AbstractClassFunction)
Compute matrix projection associated to the abstract class function `χ` and
permutation action on the set `1:d`.

The dimension `d` of the projection is equal to the degree of the permutations
in `conjugacy_classes(χ)`. Returned tuple consist of coefficient (weight) and
the matrix realization of the projection.
"""
function matrix_projection(χ::Character{T}) where {T}
    U = matrix_projection(values(χ), conjugacy_classes(χ))
    deg = degree(χ)
    ordG = sum(length, conjugacy_classes(χ))

    return U, T(deg) / ordG
end

function matrix_projection(χ::AbstractClassFunction)
    deg = degree(χ)
    if iszero(deg) # short circuting the trivial case
        return zeros(eltype(χ), 0, degree(first(first(conjugacy_classes(χ))))), deg
    end
    ordG = sum(length, conjugacy_classes(χ))
    U = matrix_projection(values(χ), conjugacy_classes(χ))

    return U, deg / ordG
end

function add_inverse_permutation!(result, val, i::Int, j::Int)
    result[i, j] += val
end

"""
    matrix_projection(vals, ccls)
Return matrix projection defined by an abstract class function with values `vals`
attained on (the corresponding) conjugacy classes `ccls`.

The dimension of the projection is equal to the degree of the permutations in `ccls`.
"""
function matrix_projection(
    vals::AbstractVector{T},
    ccls::AbstractVector{<:AbstractOrbit{<:AbstractPerm}},
    dim = degree(first(first(ccls))),
) where {T}

    result = zeros(T, dim, dim)

    for (val, cc) in zip(vals, ccls)
        for g in cc
            for i in 1:dim
                add_inverse_permutation!(result, val, i, i^g)
            end
        end
    end

    return result
end

"""
    action_character(conjugacy_cls)
Return the character of the representation given by the elements in the conjugacy classes
`conjugacy_cls`.

This corresponds to the classical definition of character as a trace of corresponding matrices.
If the action is given by permutaion, this will be an `Int`-valued Character.
"""
function action_character(conjugacy_cls::AbstractVector{<:AbstractOrbit{<:AbstractPerm}})
    vals = Int[PermutationGroups.nfixedpoints(first(cc)) for cc in conjugacy_cls]
    # in general:
    # vals = [tr(matrix_representative(first(cc))) for cc in conjugacy_cls]
    return SymbolicWedderburn.Character(vals, conjugacy_cls)
end

"""
    affordable_real!(chars::AbstractVector{<:AbstractClassFunction})
Return _real_ characters formed from `chars` by replacing `χ` with `2re(χ)` when necessary.
"""
function affordable_real!(chars::AbstractVector{T}) where {T<:AbstractClassFunction}
    pmap = PowerMap(conjugacy_classes(first(chars)))
    res = affordable_real!.(chars, Ref(pmap))
    return unique!(res)
end

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
        row_echelon_form(u_)
    else
        u .*= weight
        row_echelon_form(u)
    end
    dim = length(pivots)
    image_basis = image[1:dim, :]

    @debug "isotypical subspace corresponding to χ has dimension $dim" χ

    return image_basis
end

"""
    symmetry_adapted_basis([T::Type=Rational{Int},] G::AbstractPermutationGroup)
Compute a basis for the linear space `ℝⁿ` which is invariant under the symmetry of `G`.

The permutation group is acting naturally on `1:degree(G)`. The coefficients of
the invariant basis are returned in (orthogonal) blocks corresponding to irreducible
characters of `G`.

!!! Note:
Each block is invariant under the action of `G`, i.e. the action may permute
vectors from symmetry adapted basis within each block.
"""
symmetry_adapted_basis(G::AbstractPermutationGroup) = symmetry_adapted_basis(Rational{Int}, G)

symmetry_adapted_basis(::Type{T}, G::AbstractPermutationGroup) where {T} =
    symmetry_adapted_basis(T, characters_dixon(real(T), G))

"""
    symmetry_adapted_basis([T::Type=Rational{Int},] G::Group, basis, action)
Compute a basis for the linear space spanned by `basis` which is invariant under
the symmetry of `G`.

* The action used in these computations is
> `(b,g) → action(b,g)` for `b ∈ basis`, `g ∈ G`
and needs to be defined by the user.
* It is assumed that `G` acts on a subset of basis and the action needs to be
extended to the whole `basis`. If `G` already acts on the whole `basis`,
a call to `symmetry_adapted_basis(G)` is preferred.
* For inducing the action `basis` needs to be indexable and iterable
(e.g. in the form of an `AbstractVector`).

The coefficients of the invariant basis are returned in (orthogonal) blocks
corresponding to irreducible characters of `G`.
The `eltype` of each of those blocks will be `Cyclotomic{real(T)}`, if possible.

!!! Note:
Each block is invariant under the action of `G`, i.e. the action may permute
vectors from symmetry adapted basis within each block.
"""
symmetry_adapted_basis(G::Group, basis, action) =
    symmetry_adapted_basis(Rational{Int}, G, basis, action)

function symmetry_adapted_basis(::Type{T}, G::Group, basis, action) where {T}
    chars = characters_dixon(real(T), G)
    ehom = ExtensionHomomorphism(basis, action)
    chars_ext = let chars = chars, ehom = ehom
        ψ = ehom(first(chars))
        ext_ccG = conjugacy_classes(ψ)
        [Character(values(χ), χ.inv_of, ext_ccG) for χ in chars]
    end
    return symmetry_adapted_basis(T, chars_ext)
end

function symmetry_adapted_basis(
    ::Type{T},
    chars::AbstractVector{<:AbstractClassFunction},
) where {T}
    ψ = action_character(conjugacy_classes(first(chars)))
    real_chars = T <: Real ? affordable_real!(deepcopy(chars)) : chars

    multiplicities = Int[Int(dot(ψ, χ)) / Int(dot(χ, χ)) for χ in real_chars]
    degrees = degree.(real_chars)

    @debug "Decomposition into character spaces:
    degrees=$(join([lpad(d, 5) for d in degrees], ""))
    multips=$(join([lpad(m, 5) for m in multiplicities], ""))"

    dot(multiplicities, degrees) == degree(ψ) ||
        @warn "chars do not constitute a complete basis for action:
        $(dot(multiplicities, degrees)) ≠ $(degree(ψ))"
    constituents = [χ for (χ, m) in zip(real_chars, multiplicities) if m ≠ 0]

    return isotypical_basis.(constituents)
end
