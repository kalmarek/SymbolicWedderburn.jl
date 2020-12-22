"""
    matrix_projection(χ::AbstractClassFunction)
Compute matrix projection associated to the abstract class function `χ` and
permutation action on the set `1:d`.

The dimension `d` of the projection is equal to the degree of the permutations
in `conjugacy_classes(χ)`. Returned tuple consist of coefficient (weight) and
the matrix realization of the projection.
"""
function matrix_projection(χ::AbstractClassFunction)
    U = matrix_projection(values(χ), conjugacy_classes(χ))
    deg = degree(χ)
    ordG = sum(length, conjugacy_classes(χ))
    return deg // ordG, U
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
) where {T}

    dim = degree(first(first(ccls)))

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
    conjugate_pairs
end

"""
    real_vchars(chars::AbstractVector{<:ClassFunction})
Return _real_ virtual characters formed from `chars` by pairing conjugates with each other.

I.e. for every pair `(a, ā)` of conjugate functions is replaced by a new pair
`(a+ā, -im*(a-ā))` of real (virtual) characters.
"""
function real_vchars(chars::AbstractVector{<:AbstractClassFunction})
    pairs = _conjugate_pairs(chars)
    res = [VirtualCharacter(chars[i]) for (i, j) in pairs if i == j] # real ones
    for (i, j) in pairs
        i == j && continue
        χ, χ_bar = chars[i], chars[j]
        push!(res, χ + χ_bar)
        push!(res, -im * (χ - χ_bar))
    end

    return res
end

"""
    isotypical_basis(χ::Union{<:Character, VirtualCharacter<:})
Compute a basis of the image of the projection corresponding to a (virtual) character `χ`.

Return the coefficients of basis vectors in an invariant subspace corresponding to `χ`
(so called _isotypical subspace_) in the permutation action on `ℝⁿ` encoded by
the conjugacy classes of `χ`.
"""
function isotypical_basis(χ::ClassFunction)
    isreal(χ) || throw("Isotypical projection for complex characters is not implemented.")

    c, u = matrix_projection(χ)
    if iszero(c)
        u = similar(u, 0, size(u, 2))
    end

    # we stop doing things symbolicaly in this call to float:
    image, pivots = row_echelon_form(float.(u))
    dim = length(pivots)
    image_basis = image[1:dim, :]

    @debug "isotypical subspace corresponding to χ has dimension $dim" χ

    return image_basis
end

"""
    symmetry_adapted_basis_float(G::PermGroup)
Compute a basis for the linear space `ℝⁿ` which is invariant under the symmetry of `G`.

The group is considered to act by permutations on set `1:n`. The coefficients of
the invariant basis are returned in (orthogonal) blocks corresponding to irreducible
characters of `G`.

!!! Note:
Each block is invariant under the action of `G`, i.e. the action may permute
vectors from symmetry adapted basis within each block.
"""
symmetry_adapted_basis_float(G::PermutationGroups.PermGroup) =
    _symmetry_adapted_basis_float(characters_dixon(G))

"""
    symmetry_adapted_basis_float(G::PermGroup, basis)
Compute a basis for the linear space spanned by `basis` which is invariant under
the symmetry of `G`.

* The action used in these computations is `(b,g) → b^g` (for an element `b ∈ basis` and a group element `g ∈ G`) and needs to be defined by the user.
* It is assumed that `G` is a permutation group acting on a subset of basis and the action needs to be induced to an the whole `basis`. If `G` already acts on the whole `basis`, a call to `symmetry_adapted_basis_float(G)` is preferred.
* For inducing the action `basis` needs to be indexable and iterable (e.g. in the form of an `AbstractVector`).

The coefficients of the invariant basis are returned in (orthogonal) blocks
corresponding to irreducible characters of `G`.

!!! Note:
Each block is invariant under the action of `G`, i.e. the action may permute
vectors from symmetry adapted basis within each block.
"""
function symmetry_adapted_basis_float(G::PermutationGroups.PermGroup, basis)
    chars = characters_dixon(G)
    induced_chars = ihom.(chars)

    @debug "Double-checking the induced action..." let ccls =
            conjugacy_classes(first(chars)),
        large_gens = ihom.(gens(G))

        G_large = PermGroup(large_gens)
        ccG_large = conjugacy_classes(G_large)
        @assert all(Set.(collect.(ccG_large)) .== Set.(collect.(ccls)))
    end

    return _symmetry_adapted_basis_float(induced_chars)
end

function _symmetry_adapted_basis_float(chars::AbstractVector{<:AbstractClassFunction})
    vr_chars = real_vchars(chars)
    # return ony the non-zero blocks:
    return filter!(!iszero ∘ first ∘ size, isotypical_basis.(vr_chars))
end
