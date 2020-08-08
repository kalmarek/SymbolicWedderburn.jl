function row_echelon_form!(A::AbstractMatrix{T}) where T <: FiniteFields.GF
    c, d = size(A)
    l = Int[]
    pos = 0
    i = 0
    for i = 1:d
        j = findfirst(x -> !iszero(x), @view A[pos+1:end, i])
        j === nothing && continue
        j += pos
        pos += 1
        push!(l, i)

        if (pos != j)
            A[[pos, j], :] = A[[j, pos], :]
        end
        A[pos, :] ./= A[pos, i]
        
        for j = 1:c 
            if j!=pos
                A[j, :] -= A[j, i].*A[pos, :]
            end
        end
    end
    return A, l
end

function row_echelon_form(A::AbstractMatrix)
    return row_echelon_form!(deepcopy(A))
end

function right_nullspace(M::AbstractMatrix{T}) where T
    A, l = row_echelon_form(M)
    c, d = size(A)
    (length(l) == d) && return zeros(T, d)
    W = zeros(T, d, d-length(l))
    i = 0
    for el in setdiff(1:d, l)
        i+= 1
        W[el, i] += 1 
        for (j, k) in enumerate(l)
            if j < el
                W[k, i] -= A[j, el]
            end
        end
    end
    return W
end

function left_nullspace(M::AbstractMatrix) 
    return Matrix(transpose(right_nullspace(Matrix(transpose(M)))))
end

function left_eigen(M::AbstractMatrix{T}) where T <: FiniteFields.GF
    @assert ==(size(M)...)
    Id = Matrix{T}(I, size(M)...)
    eigen = Dict{T, typeof(M)}()
    cumdim = 0
    for i in T # brute force; TODO: use factorization of characteristic/minimal polynomials
        cumdim >= size(M, 1) && break
        #do left eigenspaces!
        basis = first(row_echelon_form!(left_nullspace(M - i*Id)))
        nullity = size(basis, 1)
        if (nullity == 1) && all(iszero, basis)
            nullity -= 1 # left_nullspace returns trivial kernel if kernel is empty
        end
        if nullity > 0
            cumdim += nullity
            eigen[i] = basis
        end
    end
    return eigen
end

function _find_l(M::AbstractMatrix)
    # this function should be redundant when defining a better structure for echelonized subspaces
    l = Int[]
    for i = 1:size(M, 2)
        j = findfirst(isone, @view M[length(l)+1:end,i])
        if j !== nothing
            push!(l, i)
        end
    end
    return l
end

# EigenSpaceDecomposition

function eigen_decomposition!(M::Matrix{T}) where T <: FiniteFields.GF
    eigspace_ptrs = Vector{Int}()
    eigen = left_eigen(M)
    sizehint!(eigspace_ptrs, length(eigen) + 1)
    push!(eigspace_ptrs, 1)
    for val in sort!(collect(keys(eigen))) #to get deterministic behaviour
        basis = eigen[val]
        dim = size(basis, 1)
        cd = eigspace_ptrs[end]
        ran = cd:cd+dim-1
        M[ran, :] = basis 
        push!(eigspace_ptrs, cd+dim)
    end
    @assert eigspace_ptrs[end] == size(M, 1) + 1
    return M, eigspace_ptrs
end

mutable struct EigenSpaceDecomposition{T <: FiniteFields.GF}
    basis::Matrix{T}
    eigspace_ptrs::Vector{Int}

    function EigenSpaceDecomposition(
                                     basis::Matrix{T},
                                     eigspace_ptrs::AbstractVector{<:Integer}
                                    ) where T <: FiniteFields.GF

        @assert eigspace_ptrs[1] == 1
        @assert eigspace_ptrs[end] == size(basis, 1) + 1
        return new{T}(basis, eigspace_ptrs)
    end
end

EigenSpaceDecomposition(M::Matrix{T}) where T <: FiniteFields.GF =
    EigenSpaceDecomposition(eigen_decomposition!(deepcopy(M))...)

function Base.show(io::IO, ::MIME"text/plain", esd::EigenSpaceDecomposition{T}) where T
    println(io, tuple(diff(esd.eigspace_ptrs)...), "-splitting over ", T)
    print(io, esd.basis)
end

function Base.show(io::IO, esd::EigenSpaceDecomposition{T}) where T
    print(io, tuple(diff(esd.eigspace_ptrs)...), "-splitting over ", T)
end

Base.length(esd::EigenSpaceDecomposition) = length(esd.eigspace_ptrs)-1

function Base.getindex(esd::EigenSpaceDecomposition, i::Int)
    @boundscheck 1 <= i <= length(esd)
    return esd.basis[esd.eigspace_ptrs[i]:esd.eigspace_ptrs[i+1]-1, :]
end

function Base.iterate(esd::EigenSpaceDecomposition, s=1)
    s > length(esd) && return nothing 
    first_last = esd.eigspace_ptrs[s]:esd.eigspace_ptrs[s+1]-1
    return (esd.basis[first_last, :], s+1)
end

LinearAlgebra.isdiag(esd::EigenSpaceDecomposition) = esd.eigspace_ptrs == 1:length(esd)+1

function refine(esd::EigenSpaceDecomposition{T}, M::Matrix{T}) where T
    nbasis = Array{T}(undef, 0, size(first(esd), 2))
    nptrs = [1]
    for (i, e) in enumerate(esd)
        if size(e, 1) > 1
            esd2, ptrs = eigen_decomposition!(e*M[:, _find_l(e)])
            nbasis = vcat(nbasis, esd2*e)
            append!(nptrs, ptrs.+(pop!(nptrs)-1))
        else
            nbasis = vcat(nbasis, e)
            push!(nptrs, nptrs[end]+1)
        end
    end
    return EigenSpaceDecomposition(nbasis,nptrs)
end
