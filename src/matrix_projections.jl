## preallocation

function _preallocate_spmatrix(::Type{T}, sizes::Tuple, sizehint) where {T}
    res = spzeros(T, sizes...)
    @static if VERSION >= v"1.7"
        sizehint!(res, first(sizes) * sizehint)
    end

    return res
end

__an_elt(χ::Character) = first(first(conjugacy_classes(χ)))
__an_elt(α::AlgebraElement) = first(basis(parent(α)))

_projection_size(p::PermutationGroups.AbstractPerm) = (d = degree(p); (d, d))
_projection_size(m::AbstractMatrix) = size(m)

_projection_size(::Nothing, χ::Character) = _projection_size(__an_elt(χ))
_projection_size(::Nothing, α::AlgebraElement) = _projection_size(__an_elt(α))
_projection_size(hom::InducedActionHomomorphism, χ::Character) =
    _projection_size(induce(hom, __an_elt(χ)))
_projection_size(hom::InducedActionHomomorphism, α::AlgebraElement) =
    _projection_size(induce(hom, __an_elt(α)))

_hint(χ::Character) = length(conjugacy_classes(χ))
_hint(α::AlgebraElement) = count(!iszero, StarAlgebras.coeffs(α))

preallocate(::Type{T}, χ::Union{Character,AlgebraElement}) where {T} =
    preallocate(T, nothing, χ)

function preallocate(
    hom::InducedActionHomomorphism,
    χ::Union{Character,AlgebraElement},
)
    T = Base._return_type(*, Tuple{coeff_type(hom), eltype(χ)})
    return preallocate(T, hom, χ)
end

function preallocate(
    ::Type{T},
    a::Union{Nothing,InducedActionHomomorphism},
    χ::Union{Character,AlgebraElement},
) where {T}
    sizes = _projection_size(a, χ)
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
matrix_projection(χ::Character{T}) where {T} =
    _mproj_outsT!(preallocate(T, χ), χ)

function matrix_projection(
    hom::InducedActionHomomorphism,
    χ::Character{T},
) where {T}
    return _mproj_outsT!(preallocate(hom, χ), hom, χ)
end

function matrix_projection(χ::Character{T}) where {T<:Rational}
    if all(isone ∘ Cyclotomics.conductor, table(χ))
        return _mproj_fitsT!(preallocate(T, χ), χ)
    else
        return _mproj_outsT!(preallocate(T, χ), χ)
    end
end

matrix_projection(χ::Character{T}) where {T<:Union{Cyclotomic,Complex}} =
    _mproj_fitsT!(preallocate(T, χ), χ)
function matrix_projection(
    hom::InducedActionHomomorphism,
    χ::Character{T},
) where {T<:Union{Cyclotomic,Complex}}
    return _mproj_fitsT!(preallocate(hom, χ), hom, χ)
end

function matrix_projection(χ::Character{T}) where {T<:AbstractFloat}
    if all(isreal, table(χ))
        return _mproj_fitsT!(preallocate(T, χ), χ)
    else
        return _mproj_outsT!(preallocate(T, χ), χ)
    end
end

function matrix_projection(
    hom::InducedActionHomomorphism,
    χ::Character{T},
) where {T<:AbstractFloat}
    if all(isreal, table(χ))
        return _mproj_fitsT!(preallocate(hom, χ), hom, χ)
    else
        return _mproj_outsT!(preallocate(hom, χ), hom, χ)
    end
end

_mproj_fitsT!(args...) = _mproj_outsT!(args...)

function _mproj_outsT!(mproj::AbstractMatrix{T}, χ::Character) where {T}
    mproj .= sum(
        c .* matrix_projection_irr(Character(table(χ), i)) for
        (i, c) in pairs(constituents(χ)) if !iszero(c)
    )
    return mproj
end

function _mproj_outsT!(
    mproj::AbstractMatrix{T},
    hom::InducedActionHomomorphism,
    χ::Character,
) where {T}
    mproj .= sum(
        c .* matrix_projection_irr(hom, Character(table(χ), i)) for
        (i, c) in pairs(constituents(χ)) if !iszero(c)
    )
    return mproj
end

"""
    matrix_projection_irr([::Type, ]χ::Character)
    matrix_projection_irr(hom::InducedActionHomomorphism, χ::Character)
Compute matrix projection associated to an *irreducible* character `χ`.

The returned matrix defines so called *isotypical* projection.

See also [matrix_projection](@ref).
"""
matrix_projection_irr(χ::Character{T}) where {T} = matrix_projection_irr(T, χ)
matrix_projection_irr(::Type{T}, χ::Character) where {T} =
    matrix_projection_irr_acc!(preallocate(T, χ), χ, 1)

matrix_projection_irr(hom::InducedActionHomomorphism, χ::Character) =
    matrix_projection_irr_acc!(preallocate(hom, χ), hom, χ, 1)

function matrix_projection_irr_acc!(
    result::AbstractMatrix,
    χ::Character,
    weight,
)
    @assert isirreducible(χ)
    LinearAlgebra.checksquare(result)
    ccls = conjugacy_classes(χ)
    vals = collect(values(χ))
    weight *= degree(χ) // sum(length, ccls)
    matrix_projection_irr_acc!(result, vals, ccls, weight)
    return result
end

function matrix_projection_irr_acc!(
    result::AbstractMatrix,
    vals,
    ccls::AbstractVector{<:AbstractOrbit{<:PermutationGroups.AbstractPerm}},
    weight,
)
    iszero(weight) && return result
    for (val, cc) in zip(vals, ccls)
        iszero(val) && continue
        w = weight * val
        for g in cc
            for i in 1:size(result, 1)
                result[i, i^g] += w
            end
        end
    end

    return result
end

function matrix_projection_irr_acc!(
    result::AbstractMatrix,
    vals,
    ccls::AbstractVector{<:AbstractOrbit{<:AbstractMatrix}},
    weight,
)
    # TODO: call to inv(Matrix(g)) is a dirty hack, since if `g`
    # is given by a sparse matrix `inv(g)` will fail.
    for (val, cc) in zip(vals, ccls)
        iszero(val) && continue
        for g in cc
            res .+= (weight * val) .* inv(convert(Matrix, g))
        end
    end
    return res
end

function matrix_projection_irr_acc!(
    result::AbstractMatrix,
    hom::InducedActionHomomorphism,
    χ::Character,
    weight,
)
    @assert isirreducible(χ)
    LinearAlgebra.checksquare(result)
    vals = collect(values(χ))
    ccls = conjugacy_classes(χ)
    weight *= degree(χ) // sum(length, ccls)
    matrix_projection_irr_acc!(result, hom, vals, ccls, weight)
    return result
end

function matrix_projection_irr_acc!(
    result::AbstractMatrix,
    hom::InducedActionHomomorphism{<:ByPermutations},
    class_values,
    conjugacy_cls,
    weight,
)
    iszero(weight) && return result
    for (val, ccl) in zip(class_values, conjugacy_cls)
        iszero(val) && continue
        w = weight * val
        for g in ccl
            h = induce(hom, g)
            @assert h isa PermutationGroups.Perm
            for i in 1:size(result, 1)
                result[i, i^h] += w
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
    weight,
)
    iszero(weight) && return result
    for (val, cc) in zip(class_values, conjugacy_cls)
        iszero(val) && continue
        w = weight * val
        for g in cc
            result .+= w .* induce(hom, inv(g))
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
    matrix_representation_acc!(preallocate(hom, α), hom, α)

function matrix_representation_acc!(
    result::AbstractMatrix,
    α::AlgebraElement{
        <:StarAlgebra{<:PermutationGroups.AbstractPermutationGroup},
    },
)
    b = basis(parent(α))
    for (idx, val) in StarAlgebras._nzpairs(StarAlgebras.coeffs(α))
        g = b[idx]
        iszero(val) && continue
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
        iszero(val) && continue
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
        iszero(val) && continue
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
* By default (symbolic, exact) row echelon form is produced, and therefore there
is no guarantee on the orthogonality of the returned basis vectors (rows).
* If `eltype(A) <: AbstractFloat` this is followed by a call to `svd` (or
`qr` in the sparse case) and the appropriate rows of its orthogonal factor are
returned, thus the basis is (numerically) orthonormal.

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
function image_basis!(A::AbstractSparseMatrix)
    A, p = row_echelon_form!(A)
    dropzeros!(A)
    return A, p
end

_eps(::Type{T}) where {T} = T <: Complex ? 2eps(real(T)) : eps(T)

function _orth!(M::AbstractMatrix{T}) where {T<:Union{AbstractFloat,Complex}}
    F = svd!(convert(Matrix{T}, M))
    M_rank = count(>(maximum(size(M)) * _eps(T)), F.S)
    return F.Vt, M_rank
end

function _orth!(
    M::AbstractSparseMatrix{T},
) where {T<:Union{AbstractFloat,Complex}}
    F = qr(M)
    M_rank = rank(F)
    result = let tmp = F.Q * Matrix(I, size(F.Q, 2), M_rank)
        sparse(tmp[invperm(F.prow), :]')
    end
    result = droptol!(result, maximum(size(result)) * _eps(T))
    return result, M_rank
end

function image_basis!(
    A::AbstractMatrix{T},
) where {T<:Union{AbstractFloat,Complex}}
    A, p = row_echelon_form!(A)
    A_orth, A_rank = _orth!(@view A[1:length(p), :])
    @assert A_rank == length(p) "_orth rank doesn't agree with rref rank!"
    return A_orth, 1:A_rank
end

function image_basis!(
    A::AbstractSparseMatrix{T},
) where {T<:Union{AbstractFloat,Complex}}
    N = LinearAlgebra.checksquare(A)
    droptol!(A, N * _eps(T))
    # A, p = row_echelon_form!(A)
    A_orth, A_rank = _orth!(A)
    # @assert A_rank == length(p) "_orth rank doesn't agree with rref rank!"
    return A_orth, 1:A_rank
end
