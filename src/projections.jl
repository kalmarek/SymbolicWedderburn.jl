function matrix_projection(χ::AbstractClassFunction)
    U = matrix_projection(values(χ), conjugacy_classes(χ))
    deg = degree(χ)
    ordG = sum(length, conjugacy_classes(χ))
    return deg // ordG, U
end

function matrix_projection(
    vals::AbstractVector{T},
    ccls::AbstractVector{<:AbstractOrbit{<:Perm}},
    d::Integer = degree(first(first(ccls))),
) where {T}

    result = zeros(T, d, d)

    for (val, cc) in zip(vals, ccls)
        for g in cc
            for i = 1:d
                result[i, i^g] += val
            end
        end
    end

    return result
end

function _conjugate_pairs(chars::Vector{<:AbstractClassFunction})
    vals = values.(chars)
    tovisit = trues(length(vals))
    conjugate_pairs = Dict{Int, Int}()
    for (i, v) in enumerate(vals)
        tovisit[i] || continue
        k = findfirst(==(conj(v)), vals)
        conjugate_pairs[i] = k
        tovisit[i] = false
        tovisit[k] = false
    end
    conjugate_pairs
end

function real_vchars(chars::AbstractVector{<:AbstractClassFunction})
    pairs = _conjugate_pairs(chars)
    res = [VirtualCharacter(chars[i]) for (i,j) in pairs if i == j] # real ones
    for (i,j) in pairs
        i == j && continue
        χ, χ_bar = chars[i], chars[j]
        push!(res, χ + χ_bar)
        push!(res, -E(4)*(χ - χ_bar))
    end

    return res
end

"""
    isotypical_basis(χ::Union{<:Character, VirtualCharacter<:})
Compute a basis of the image of the projection corresponding to a (virtual) character `χ`.

Return the coefficients of basis vectors in an invariant subspace corresponding to `χ`
(so called _isotypical subspace_) in the permutation action on `ℝⁿ`
encoded by the conjugacy classes of `χ`.
"""
function isotypical_basis(χ::ClassFunction)
    isreal(χ) || throw("Isotypical projection for complex characters is not implemented.")

    c, u = matrix_projection(χ)
    if iszero(c)
        u = similar(u, 0, size(u, 2))
    end

    # we stop doing things symbolicaly in this call to float:
    image, pivots = SymbolicWedderburn.row_echelon_form(float.(u))
    dim = length(pivots)
    image_basis = image[1:dim, :]

    @debug "isotypical subspace corresponding to χ has dimension $dim" χ

    return image_basis
end

"""
    symmetry_adapted_basis(G::Group, basis)
Compute a basis for the linear space spanned by `basis` which is invariant under
the symmetry of `G`. The action used in these computations is `(b,g) → b^g`
(for an element `b ∈ basis` and a group element `g ∈ G`) and needs to be defined by the user.

Return coefficients of the invariant basis in (orthogonal) blocks corresponding to irreducible characters of `G`.

!!! Note:
Each block is invariant under the action of `G`, i.e. the action may permute vectors from
symmetry adapted basis within each block.
"""
function symmetry_adapted_basis(G::PermutationGroups.Group, basis)
    chars = characters_dixon(G)

    @debug "Characters share conjugacy classes" @assert all(
        χ.inv_of == first(chars).inv_of for χ in chars
    )

    h = InducingHomomorphism(basis)
    ccG_large = h.(conjugacy_classes(first(chars)))

    @debug "Double-checking the induced action..." let ccls = ccG_large,
        large_gens = h.(gens(G))

        G_large = PermGroup(large_gens)
        ccG_large = conjugacy_classes(G_large)
        @assert all(Set.(collect.(ccG_large)) .== Set.(collect.(ccls)))
    end

    chars_action_on_basis = [Character(values(χ), χ.inv_of, ccG_large) for χ in chars]
    vr_chars = real_vchars(chars_action_on_basis)

    # return ony the non-zero subbases:
    return filter!(!iszero ∘ first ∘ size, isotypical_basis.(vr_chars))
end
