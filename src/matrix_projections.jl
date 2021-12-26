## preallocation

function _preallocate_spmatrix(::Type{T}, sizes, sizehint) where {T}
    res = spzeros(T, sizes)
    sizehint!(res, first(sizes) * sizehint)
    return res
end

__an_elt(χ::Character) = first(first(conjugacy_classes(χ)))
__an_elt(α::AlgebraElement) = first(basis(parent(α)))

_projection_size(p::PermutationGroups.AbstractPerm) = (d = degree(p); (d, d))
_projection_size(m::AbstractMatrix) = size(m)

_projection_size(χ::Character) = _projection_size(__an_elt(χ))
_projection_size(hom::InducedActionHomomorphism, χ::Character) =
    _projection_size(induce(hom, __an_elt(χ)))
_projection_size(α::AlgebraElement) = _projection_size(__an_elt(α))
_projection_size(hom::InducedActionHomomorphism, α::AlgebraElement) =
    _projection_size(induce(hom, __an_elt(χ)))

_hint(χ::Character) = length(conjugacy_classes(χ))
_hint(α::AlgebraElement) = count(!iszero, StarAlgebras.coeffs(α))

function preallocate(::Type{T}, χ::Union{Character, AlgebraElement}) where {T}
    sizes = _projection_size(χ)
    return _preallocate_spmatrix(T, sizes, _hint(χ))
end

function preallocate(
    ::Type{T},
    hom::InducedActionHomomorphism,
    χ::Union{Character, AlgebraElement},
) where {T}
    sizes = _projection_size(hom, χ)
    return _preallocate_spmatrix(T, sizes, _hint(χ))
end

## matrix projection [irreducible]
"""
    matrix_projection([hom::InducedActionHomomorphism, ]χ::Character{T})
Compute matrix projection associated to character `χ`.

Returned `M<:AbstractMatrix{T}` of size `(d, d)` where the degree `d` of the
projecion is determined by elements in `conjugacy_classes(χ)`. E.g. `d` could
be equal to the `degree` when conjugacy classes consist of `AbstractPerms`.
If the homomorphism is passed, the dimension will be derived in similar
manner from the elements of the image of the homomorphism.

The precise type of `M` can be altered by overloading

```
preallocate(::Type{T}, [hom::InducedActionHomomorphism, ]χ::Character)
```
"""
function matrix_projection(χ::Character{T}) where {T}
    tbl = table(χ)
    mproj = preallocate(T, χ)
    for (c, ψ) in zip(constituents(χ), irreducible_characters(T, tbl))
        iszero(c) && continue
        mproj .= c .* matrix_projection_irr_acc!(mproj, ψ)
    end
    return mproj
end

function matrix_projection(
    hom::InducedActionHomomorphism,
    χ::Character{T},
) where {T}
    tbl = table(χ)
    mproj = preallocate(T, hom, χ)

    for (c, ψ) in zip(constituents(χ), irreducible_characters(T, tbl))
        iszero(c) && continue
        mproj .= c .* matrix_projection_irr_acc!(mproj, hom, ψ)
    end
    return mproj
end

"""
    matrix_projection_irr(χ::Character)
    matrix_projection_irr([::Type, hom::InducedActionHomomorphism, ]χ::Character)
Compute matrix projection associated to an *irreducible* character `χ`.

The returned matrix defines so called *isotypical* projection.

See also [matrix_projection](@ref).
"""
matrix_projection_irr(χ::Character{T}) where {T} = matrix_projection_irr(T, χ)

matrix_projection_irr(::Type{T}, χ::Character) where {T} =
    matrix_projection_irr_acc!(preallocate(T, χ), χ)

function matrix_projection_irr(
    hom::InducedActionHomomorphism,
    χ::Character{T},
) where {T}
    return matrix_projection_irr(T, hom, χ)
end

function matrix_projection_irr(
    ::Type{T},
    hom::InducedActionHomomorphism,
    χ::Character,
) where {T}
    return matrix_projection_irr_acc!(preallocate(T, hom, χ), hom, χ)
end

function matrix_projection_irr_acc!(result::AbstractMatrix, χ::Character)
    @assert isirreducible(χ)
    LinearAlgebra.checksquare(result)
    ccls = conjugacy_classes(χ)
    vals = collect(values(χ))
    result .= matrix_projection_irr_acc!(result, vals, ccls)
    result .*= degree(χ) // sum(length, ccls)
    return result
end

function matrix_projection_irr_acc!(
    result::AbstractMatrix,
    vals,
    ccls::AbstractVector{<:AbstractOrbit{<:PermutationGroups.AbstractPerm}},
) where {T}
    for (val, cc) in zip(vals, ccls)
        for g in cc
            for i in 1:size(result, 1)
                result[i, i^g] += val
            end
        end
    end

    return result
end

function matrix_projection_irr_acc!(
    result::AbstractMatrix,
    vals,
    ccls::AbstractVector{<:AbstractOrbit{<:AbstractMatrix}},
)
    # TODO: call to inv(Matrix(g)) is a dirty hack, since if `g`
    # is given by a sparse matrix `inv(g)` will fail.
    for (val, cc) in zip(vals, ccls)
        for g in cc
            res .+= val .* inv(convert(Matrix, g))
        end
    end
    return res
end

function matrix_projection_irr_acc!(
    result::AbstractMatrix,
    hom::InducedActionHomomorphism,
    χ::Character,
)
    @assert isirreducible(χ)
    LinearAlgebra.checksquare(result)
    ccls = conjugacy_classes(χ)
    vals = collect(values(χ))
    result .= matrix_projection_irr_acc!(result, hom, vals, ccls)
    result .*= degree(χ) // sum(length, ccls)
    return result
end

function matrix_projection_irr_acc!(
    result::AbstractMatrix,
    hom::InducedActionHomomorphism{<:ByPermutations},
    class_values,
    conjugacy_cls,
)
    for (val, ccl) in zip(class_values, conjugacy_cls)
        for g in ccl
            h = induce(hom, g)
            @assert h isa PermutationGroups.Perm
            for i in 1:size(result, 1)
                result[i, i^h] += val
            end
        end
    end
    return result
end

function matrix_projection_irr_acc!(
    result::AbstractMatrix,
    hom::InducedActionHomomorphism{<:ByLinearTransformation},
    class_values,
    conjugacy_cls,
)
    for (val, cc) in zip(class_values, conjugacy_cls)
        iszero(val) && continue
        for g in cc
            result .+= val .* induce(hom, inv(g))
        end
    end
    return result
end

"""
    matrix_representation([::Type{T}, ]α::AlgebraElement)
    matrix_representation([::Type{T}, ]hom::InducedActionHomomorphism, α::AlgebraElement)
Compute matrix representative of the given algebra element `α`.

If `α` is an element of the group (or monoid) algebra `RG` of `G`, then every
linear representation `π : G → GL(V)` gives rise to an algebra homomorphism
`π : RG → B(V)` (bounded operators on `V`). This function evaluates this
representation given either by
 * the natural action of `G` on some `V` (e.g. every permutation group of degree
 `d` acts on `d`-dimensional `V` by basis vector permutation), or
 * by action homomorphism `hom`.

See also [matrix_projection](@ref).
"""
matrix_representation(α::AlgebraElement) = matrix_representation(eltype(α), α)
matrix_representation(::Type{T}, α::AlgebraElement) where {T} =
    matrix_representation_acc!(preallocate(T, α), α)
matrix_representation(hom::InducedActionHomomorphism, α::AlgebraElement) =
    matrix_representation(eltype(α), hom, α)

function matrix_representation(
    ::Type{T},
    hom::InducedActionHomomorphism,
    α::AlgebraElement,
) where {T}
    return matrix_representation_acc!(preallocate(T, hom, α), hom, α)
end

function matrix_representation_acc!(
    result::AbstractMatrix,
    α::AlgebraElement{
        <:StarAlgebra{<:PermutationGroups.AbstractPermutationGroup},
    },
)
    b = basis(parent(α))
    for (idx, val) in StarAlgebras._nzpairs(StarAlgebras.coeffs(α))
        g = b[idx]
        for i in 1:size(result, 1)
            result[i, i^g] += val
        end
    end

    return result
end

function matrix_representation_acc!(
    result::AbstractMatrix,
    hom::InducedActionHomomorphism{<:ByPermutations},
    α::AlgebraElement,
)
    b = basis(parent(α))
    for (idx, val) in StarAlgebras._nzpairs(StarAlgebras.coeffs(α))
        g = induce(hom, b[idx])
        @assert g isa PermutationGroups.Perm
        for i in 1:size(result, 1)
            result[i, i^g] += val
        end
    end

    return result
end

function matrix_representation_acc!(
    result::AbstractMatrix,
    hom::InducedActionHomomorphism{<:ByLinearTransformation},
    α::AlgebraElement,
)
    b = basis(parent(α))
    for (idx, val) in StarAlgebras._nzpairs(StarAlgebras.coeffs(α))
        result .+= val .* induce(hom, inv(b[idx]))
    end

    return result
end

## Finding basis of the row-space (right image) of an AbstractMatrix

"""
    image_basis(A::AbstractMatrix)
    image_basis([hom::InducedActionHomomorphism, ]χ::Character)
    image_basis([hom::InducedActionHomomorphism, ]α::AlgebraElement)
Return basis of the row-space of a matrix.

For characters or algebra elements return a basis of the row-space of
`matrix_projection([hom, ]χ)` or `matrix_representation([hom, ]α)`.

Two methods are employed to achieve the goal.
* By default (symbolic, exact) row echelon form is produced, and therefore there is no guarantee on the orthogonality of the returned basis vectors (rows).
* If `eltype(A) <: AbstractFloat` this is followed by a call to `svd` and the appropriate rows of its `Vt` are returned, thus the basis is (numerically) orthonormal.

# Examples:
```julia
julia> a = rand(-3:3, 3,3).//rand(1:4, 3,3);

julia> a[3, :] .= 2a[1, :] .- 1a[2, :]; a
3×3 Matrix{Rational{Int64}}:
  3//2   0//1  1//1
 -1//2  -2//1  1//2
  7//2   2//1  3//2

julia> ib = SymbolicWedderburn.image_basis(a)
2×3 Matrix{Rational{Int64}}:
 1//1  0//1   2//3
 0//1  1//1  -5//12

julia> ibf = SymbolicWedderburn.image_basis(float.(a))
2×3 Matrix{Float64}:
 -0.666651  0.416657  -0.618041
  0.529999  0.847998  -8.85356e-17

```
"""
image_basis(A::AbstractMatrix) =
    ((m, p) = image_basis!(deepcopy(A)); m[1:length(p), :])

function image_basis(χ::Character)
    mpr = isirreducible(χ) ? matrix_projection_irr(χ) : matrix_projection(χ)
    image, pivots = image_basis!(mpr)
    return image[1:length(pivots), :]
end

function image_basis(hom::InducedActionHomomorphism, χ::Character)
    mpr =
        isirreducible(χ) ? matrix_projection_irr(hom, χ) :
        matrix_projection(hom, χ)
    image, pivots = image_basis!(mpr)
    return image[1:length(pivots), :]
end

function image_basis(α::AlgebraElement)
    image, pivots = image_basis!(matrix_representation(α))
    return image[1:length(pivots), :]
end

function image_basis(hom::InducedActionHomomorphism, α::AlgebraElement)
    image, pivots = image_basis!(matrix_representation(hom, α))
    return image[1:length(pivots), :]
end

image_basis!(A::AbstractMatrix) = row_echelon_form!(A)

function image_basis!(A::AbstractMatrix{T}) where {T<:AbstractFloat}
    A, p = row_echelon_form!(A)
    fact = svd!(convert(Matrix{T}, @view(A[1:length(p), :])))
    A_rank = sum(fact.S .> maximum(size(A)) * eps(T))
    @assert A_rank == length(p)
    return fact.Vt, 1:length(p)
end

function image_basis!(A::AbstractMatrix{T}) where {T<:Complex}
    fact = svd!(A)
    A_rank = sum(fact.S .> maximum(size(A)) * 2eps(real(T)))
    return fact.Vt, 1:A_rank
end
