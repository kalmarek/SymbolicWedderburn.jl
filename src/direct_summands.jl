struct DirectSummand{T,M<:AbstractMatrix{T}} <: AbstractMatrix{T}
    basis::M
    multiplicity::Int
    degree::Int

    function DirectSummand(
        basis::M,
        multiplicity::Integer,
        degree::Integer,
    ) where {M<:AbstractMatrix}
        @assert size(basis, 1) < size(basis, 2)
        # @assert 1 ≤ degree ≤ size(basis, 1)
        let (pr_rank, r) = divrem(size(basis, 1), multiplicity)
            @assert 1 ≤ pr_rank ≤ degree
            @assert r == 0
        end
        return new{eltype(M),M}(basis, multiplicity, degree)
    end
end

image_basis(ds::DirectSummand) = ds.basis
PermutationGroups.degree(ds::DirectSummand) = ds.degree
multiplicity(ds::DirectSummand) = ds.multiplicity
pr_rank(ds::DirectSummand) = div(size(image_basis(ds), 1), multiplicity(ds))
issimple(ds::DirectSummand) = isone(pr_rank(ds))

Base.size(ds::DirectSummand) = size(ds.basis)
Base.@propagate_inbounds function Base.getindex(ds::DirectSummand, i...)
    return image_basis(ds)[i...]
end

function SparseArrays.sparse(ds::DirectSummand)
    sp = sparse(image_basis(ds))
    return DirectSummand(sp, multiplicity(ds), degree(ds), issimple(ds))
end

function SparseArrays.droptol!(
    ds::DirectSummand{T,<:AbstractSparseArray},
    tol,
) where {T}
    droptol!(image_basis(ds), tol)
    return ds
end
