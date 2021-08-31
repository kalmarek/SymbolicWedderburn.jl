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
    ccls::AbstractVector{<:AbstractOrbit{<:AbstractPerm}},
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

function matrix_projection_irr(hom::InducedActionHomomorphism, χ::Character)
    @assert isirreducible(χ)
    mproj = matrix_projection_irr(hom, collect(values(χ)), conjugacy_classes(χ))
    mproj .*= degree(χ) // sum(length, conjugacy_classes(χ)) # χ(1)/order(G)
    return mproj
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
    dim = length(features(hom)),
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
    α::AlgebraElement{<:StarAlgebra{<:PermutationGroups.AbstractPermutationGroup}},
    dim = degree(parent(parent(α)))
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
