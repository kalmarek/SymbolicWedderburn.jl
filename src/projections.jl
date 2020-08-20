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

