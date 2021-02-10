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
        return zeros(eltype(χ), 0, degree(first(first(conjugacy_classes(χ))))),
        deg
    end
    ordG = sum(length, conjugacy_classes(χ))
    U = matrix_projection(values(χ), conjugacy_classes(χ))

    return U, deg / ordG
end

"""
    matrix_projection(vals, ccls)
Return matrix projection defined by an abstract class function with values `vals`
attained on (the corresponding) conjugacy classes `ccls`.

The dimension of the projection is equal to the degree of the permutations in `ccls`.
"""
function matrix_projection(
    vals::AbstractVector{T},
    ccls::AbstractVector{<:AbstractOrbit{<:Perm}},
    dim = degree(first(first(ccls))),
) where {T}

    result = zeros(T, dim, dim)

    for (val, cc) in zip(vals, ccls)
        for g in cc
            for i = 1:dim
                result[i, i^g] += val
            end
        end
    end

    return result
end

function _conjugate_pairs(chars::AbstractVector{<:ClassFunction})
    vals = values.(chars)
    tovisit = trues(length(vals))
    conjugate_pairs = Dict{Int,Int}()
    for (i, v) in enumerate(vals)
        tovisit[i] || continue
        tovisit[i] = false
        if all(isreal, v)
            tovisit[i] = false
            conjugate_pairs[i] = i
        else
            k = findfirst(==(conj(v)), vals)
            conjugate_pairs[i] = k
            tovisit[k] = false
        end
    end
    return conjugate_pairs
end

"""
    _real_vchars(chars::AbstractVector{<:AbstractClassFunction})
Return _real_ virtual characters formed from `chars` by pairing conjugates with each other.

That is for every pair `(a, ā)` of conjugate functions is replaced by a new pair
`( (a+ā)/2, -im*(a-ā)/2 )` of real (virtual) characters.
"""
function _real_vchars(chars::AbstractVector{<:AbstractClassFunction})
    pairs = _conjugate_pairs(chars)

    res = [VirtualCharacter(chars[i]) for (i, j) in pairs if i == j] # real ones
    for (i, j) in pairs
        i == j && continue
        χ, χ_bar = chars[i], chars[j]
        push!(res, (χ + χ_bar) / 2)
        push!(res, -im * (χ - χ_bar) / 2)
    end

    return res
end

"""
    isotypical_basis(χ::AbstractClassFunction)
Compute a basis of the image of the projection corresponding to a class function `χ`.

Return the coefficients of basis vectors in an invariant subspace corresponding to `χ`
(so called _isotypical subspace_) in the action on `ℝⁿ` encoded by the conjugacy classes of `χ`.
"""
function isotypical_basis(χ::AbstractClassFunction) where {T}

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
    symmetry_adapted_basis([T::Type=Rational{Int},] G::PermGroup)
Compute a basis for the linear space `ℝⁿ` which is invariant under the symmetry of `G`.

The group is considered to act by permutations on set `1:n`. The coefficients of
the invariant basis are returned in (orthogonal) blocks corresponding to irreducible
characters of `G`.

!!! Note:
Each block is invariant under the action of `G`, i.e. the action may permute
vectors from symmetry adapted basis within each block.
"""
symmetry_adapted_basis(G::PermutationGroups.PermGroup) =
    symmetry_adapted_basis(Rational{Int}, G)
symmetry_adapted_basis(
    ::Type{T},
    G::PermutationGroups.PermGroup,
) where {T<:Real} = _real_symmetry_adapted_basis(characters_dixon(T, G))
symmetry_adapted_basis(
    ::Type{T},
    G::PermutationGroups.PermGroup,
) where {T<:Complex} =
    _complex_symmetry_adapted_basis(characters_dixon(real(T), G))

"""
    symmetry_adapted_basis([T::Type=Rational{Int},] G::PermGroup, basis)
Compute a basis for the linear space spanned by `basis` which is invariant under
the symmetry of `G`.

* The action used in these computations is `(b,g) → b^g` (for an element `b ∈ basis`
and a group element `g ∈ G`) and needs to be defined by the user.
* It is assumed that `G` acts by permutations on a subset of basis and the action
needs to be extended to the whole `basis`. If `G` already acts on the whole `basis`,
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
symmetry_adapted_basis(G::PermutationGroups.PermGroup, basis) =
    symmetry_adapted_basis(Rational{Int}, G::PermutationGroups.PermGroup, basis)

function symmetry_adapted_basis(
    ::Type{T},
    G::AbstractPermutationGroup,
    basis,
) where {T<:Number}
    chars = characters_dixon(real(T), G)
    ehom = ExtensionHomomorphism(basis)
    chars_ext = ehom.(chars)

    @debug "Double-checking the induced action..." let ccls =
            conjugacy_classes(first(chars)),
        large_gens = ehom.(gens(G))

        G_large = PermGroup(large_gens)
        ccG_large = conjugacy_classes(G_large)
        @assert all(Set.(collect.(ccG_large)) .== Set.(collect.(ccls)))
    end

    if T <: Real
        return _real_symmetry_adapted_basis(chars_ext)
    else # if T <: Complex
        return _complex_symmetry_adapted_basis(chars_ext)
    end
end

function _complex_symmetry_adapted_basis(chars)
    return filter!(
        !iszero ∘ first ∘ size,
        isotypical_basis.(VirtualCharacter.(chars)),
    )
end

function _real_symmetry_adapted_basis(chars)
    return filter!(
        !iszero ∘ first ∘ size,
        isotypical_basis.(_real_vchars(chars)),
    )
end
