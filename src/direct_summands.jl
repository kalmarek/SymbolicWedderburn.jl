struct DirectSummand{T, M<:AbstractMatrix{T}} <: AbstractMatrix{T}
    basis::M
    multiplicity::Int
    degree::Int
    simple::Bool
end

StarAlgebras.basis(ds::DirectSummand) = ds.basis
issimple(ds::DirectSummand) = ds.simple
PermutationGroups.degree(ds::DirectSummand) = ds.degree
multiplicity(ds::DirectSummand) = ds.multiplicity

Base.size(ds::DirectSummand) = size(ds.basis)
Base.@propagate_inbounds Base.getindex(ds::DirectSummand, i...) =
    basis(ds)[i...]

function SparseArrays.sparse(ds::DirectSummand)
    sp = sparse(basis(ds))
    return DirectSummand(sp, multiplicity(ds), degree(ds), issimple(ds))
end

function SparseArrays.droptol!(
    ds::DirectSummand{T, <:AbstractSparseArray},
    tol
) where T
    droptol!(basis(ds), tol)
    return ds
end
