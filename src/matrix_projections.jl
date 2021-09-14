"""
    matrix_projection_irr([hom::InducedActionHomomorphism, ]χ::Character)
Compute matrix projection associated to the irreducible character `χ`.

If the homomorphism is not passed, the dimension `d` of the projection is
derived from `conjugacy_classes(χ)`. E.g. if `conjugacy_classes(χ)` consist
of permutations of degree `d` (i.e. acting naturally on the set `1:d`) the
result will be a matrix of size `(d,d)`.

If the homomorphism is passed, the dimension will be derived in similar
manner from the elements of the image of the homomorphism.
"""
function matrix_projection_irr(χ::Character)
    @assert isirreducible(χ)
    mproj = matrix_projection_irr(collect(values(χ)), conjugacy_classes(χ))
    mproj .*= degree(χ) // sum(length, conjugacy_classes(χ)) # χ(1)/order(G)
    return mproj
end

function matrix_projection_irr(
    vals,
    ccls::AbstractVector{<:AbstractOrbit{<:PermutationGroups.AbstractPerm}},
)
    dim = degree(first(first(ccls)))
    result = zeros(eltype(vals), dim, dim)

    for (val, cc) in zip(vals, ccls)
        for g in cc
            for i in 1:dim
                result[i, i^g] += val
            end
        end
    end

    return result
end

function matrix_projection_irr(
    vals,
    ccls::AbstractVector{<:AbstractOrbit{<:AbstractMatrix}},
)
    # TODO: call to inv(Matrix(g)) is a dirty hack, since if `g`
    # is given by a sparse matrix `inv(g)` will fail.
    return sum(val .* sum(g -> inv(Matrix(g)), cc) for (val, cc) in zip(vals, ccls))
end

## versions with InducedActionHomomorphisms

function matrix_projection_irr(hom::InducedActionHomomorphism, χ::Character{T}) where T
    @assert isirreducible(χ)
    mproj = matrix_projection_irr(hom, collect(values(χ)), conjugacy_classes(χ))
    mproj .*= degree(χ) // sum(length, conjugacy_classes(χ)) # χ(1)/order(G)
    return eltype(mproj) == T ? mproj : T.(mproj)
end

function matrix_projection_irr(
    hom::InducedActionHomomorphism{<:ByPermutations},
    class_values,
    conjugacy_cls,
    )
    dim = degree(induce(hom, first(first(conjugacy_cls))))
    result = zeros(eltype(class_values), dim, dim)
    for (val, ccl) in zip(class_values, conjugacy_cls)
        for g in ccl
            h = induce(hom, g)
            for i in 1:dim
                result[i, i^h] += val
            end
        end
    end
    return result
end

function matrix_projection_irr(
    hom::InducedActionHomomorphism{<:ByLinearTransformation},
    class_values,
    conjugacy_cls,
)
    return sum(
        val .* sum(g -> induce(hom, inv(g)), cc)
        for (val, cc) in zip(class_values, conjugacy_cls)
    )
end

function matrix_projection(χ::Character{T}) where T
    tbl = table(χ)
    res = sum(
        c .* matrix_projection_irr(ψ)
        for (c, ψ) in zip(constituents(χ), irreducible_characters(tbl)) if !iszero(c)
    )
    return eltype(res) == T ? res : T.(res)
end

function matrix_projection(
    hom::InducedActionHomomorphism{<:ByPermutations},
    α::AlgebraElement,
    dim::Integer = length(features(hom)),
)
    result = zeros(eltype(α), dim, dim)
    b = basis(parent(α))

    @inbounds for (j, a) in StarAlgebras._nzpairs(StarAlgebras.coeffs(α))
        g = induce(hom, b[j])
        for i in 1:dim
            result[i, i^g] += a
        end
    end

    return result
end

function matrix_projection(
    hom::InducedActionHomomorphism{<:ByLinearTransformation},
    α::AlgebraElement,
    dim::Integer = length(features(hom)),
)

    result = zeros(Base._return_type(*, Tuple{eltype(α), coeff_type(action(hom))}), dim, dim)
    b = basis(parent(α))

    for (j, a) in StarAlgebras._nzpairs(StarAlgebras.coeffs(α))
        result += a*induce(hom, inv(b[j]))
    end

    return result
end

function matrix_projection(
    α::AlgebraElement{<:StarAlgebra{<:PermutationGroups.AbstractPermutationGroup}},
    dim::Integer = degree(parent(parent(α)))
)
    result = zeros(eltype(α), dim, dim)
    b = basis(parent(α))

    @inbounds for (j, a) in StarAlgebras._nzpairs(StarAlgebras.coeffs(α))
        g = b[j]
        for i in 1:dim
            result[i, i^g] += a
        end
    end

    return result
end

image_basis!(A::AbstractMatrix) = row_echelon_form!(A)

function image_basis!(A::AbstractMatrix{T}) where {T<:AbstractFloat}
    A, p = row_echelon_form!(A)
    fact = svd!(Matrix(@view A[1:length(p), :]))
    A_rank = sum(fact.S .> maximum(size(A)) * eps(T))
    return fact.Vt, 1:A_rank
end

function image_basis!(A::AbstractMatrix{T}) where {T<:Complex}
    fact = svd!(A)
    A_rank = sum(fact.S .> maximum(size(A)) * 2eps(real(T)))
    return fact.Vt, 1:A_rank
end

image_basis(A::AbstractMatrix) =
    ((m, p) = image_basis!(deepcopy(A)); m[1:length(p), :])

function image_basis(χ::Character)
    mpr = isirreducible(χ) ? matrix_projection_irr(χ) : matrix_projection(χ)
    image, pivots = image_basis!(mpr)
    return image[1:length(pivots), :]
end

function image_basis(hom::InducedActionHomomorphism, χ::Character)
    mpr = isirreducible(χ) ? matrix_projection_irr(hom, χ) : matrix_projection(hom, χ)
    image, pivots = image_basis!(mpr)
    return image[1:length(pivots), :]
end

function image_basis(α::AlgebraElement)
    image, pivots = image_basis!(matrix_projection(α))
    return image[1:length(pivots), :]
end

function image_basis(hom::InducedActionHomomorphism, α::AlgebraElement)
    image, pivots = image_basis!(matrix_projection(hom, α))
    return image[1:length(pivots), :]
end
