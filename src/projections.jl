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
