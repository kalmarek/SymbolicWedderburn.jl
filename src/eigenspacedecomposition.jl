# Version over exact field:
function row_echelon_form!(A::AbstractMatrix{T}) where {T<:Union{FiniteFields.GF, Cyclotomic{<:Rational}, Rational}}
    pivots = Int[]
    pos = 0
    for i = 1:size(A, 2)
        j = findnext(x -> !iszero(x), @view(A[:, i]), pos+1)
        j === nothing && continue
        pos += 1
        push!(pivots, i)

        # swap the rows so that
        if (pos != j)
            for k in 1:size(A, 2)
                A[pos, k], A[j, k] = A[j, k], A[pos, k]
            end
        end
        # A[pos, :] is the row with leading nz in i-th column

        w = inv(A[pos, i])
        # all columns to the left of i are zeroed (below pivots) so we start at i
        for colidx in i:size(A, 2)
            A[pos, colidx] *= w
        end
        # A[pos, i] is 1 now

        # zero the whole i-th column above and below pivot:
        for rowidx = 1:size(A, 1)
            rowidx == pos && continue # don't zero the active row
            v = A[rowidx, i]
            iszero(v) && continue
            for k in i:size(A, 2) # to the left of pivot everything is zero
                A[rowidx, k] -= v * A[pos, k]
            end
        end
    end
    return A, pivots
end

function row_echelon_form!(A::AbstractMatrix{C}) where {T<:AbstractFloat, C<:Cyclotomic{T}}
    pivots = Int[]
    pos = 0
    for i = 1:size(A, 2)

        j = findnext(x -> abs(x) > sqrt(eps(T)), @view(A[:, i]), pos+1)
        j === nothing && continue
        pos += 1
        push!(pivots, i)

        if (pos != j)
            for k in 1:size(A, 2)
                A[pos, k], A[j, k] = A[j, k], A[pos, k]
            end
        end

        @assert abs(A[pos, i]) >= sqrt(eps(T))
        w = inv(A[pos, i])

        for colidx in i:size(A, 2)
            if abs(A[pos, colidx]) <= sqrt(eps(T))
                Cyclotomics.zero!(A[pos, colidx])
            elseif colidx == i
                Cyclotomics.one!(A[pos, colidx])
            else
                A[pos, colidx] *= w
            end
        end

        for rowidx = 1:size(A, 1)
            rowidx == pos && continue
            v = A[rowidx, i]
            if abs(v) < sqrt(eps(T))
                Cyclotomics.zero!(v)
            else
                for k in i:size(A, 2)
                    A[rowidx, k] -= v * A[pos, k]
                end
            end
        end
    end

    for i in eachindex(A)
        ndigs = floor(Int, -log10(eps(T)) - log10(max(size(A)...)) - 3)
        A[i] = Cyclotomics.roundcoeffs!(A[i], digits=ndigs)
    end

    return A, pivots
end

row_echelon_form(A::AbstractMatrix) = row_echelon_form!(deepcopy(A))

image_basis(A::AbstractMatrix) = image_basis!(deepcopy(A))

image_basis!(A::AbstractMatrix) = row_echelon_form!(A)

function image_basis!(A::AbstractMatrix{T}) where {T<:AbstractFloat}
    fact = svd!(A)
    A_rank = sum(fact.S .> maximum(size(A)) * eps(T))
    return fact.Vt, 1:A_rank
end

function image_basis!(A::AbstractMatrix{T}) where {T<:Complex}
    fact = svd!(A)
    A_rank = sum(fact.S .> maximum(size(A)) * 2eps(real(T)))
    return fact.Vt, 1:A_rank
end

function right_nullspace(M::AbstractMatrix{T}) where {T}
    A, pivots = row_echelon_form(M)
    ncolsA = size(A, 2)
    length(pivots) == ncolsA && return zeros(T, ncolsA, 0)
    W = zeros(T, ncolsA, ncolsA - length(pivots))
    for (i, el) in enumerate(setdiff(1:ncolsA, pivots))
        W[el, i] += 1
        for (j, k) in enumerate(pivots)
            if j < el
                W[k, i] -= A[j, el]
            end
        end
    end
    return W
end

function left_nullspace(M::AbstractMatrix)
    return transpose(right_nullspace(transpose(M)))
end


function left_eigen(M::AbstractMatrix{T}) where {T<:FiniteFields.GF}
    @assert ==(size(M)...)
    eigen = Dict{T,typeof(M)}()
    cumdim = 0
    for i in T # brute force; TODO: use factorization of characteristic/minimal polynomials
        cumdim >= size(M, 1) && break
        #do left eigenspaces!
        basis = first(row_echelon_form!(left_nullspace(M - i * I)))
        if size(basis, 1) > 0 # nullity
            cumdim += size(basis, 1)
            eigen[i] = basis
        end
    end
    return eigen
end

function _find_l(M::AbstractMatrix)
    # this function should be redundant when defining a better structure for echelonized subspaces
    l = Int[]
    for i = 1:size(M, 2)
        j = findfirst(isone, @view M[length(l)+1:end, i])
        if j !== nothing
            push!(l, i)
        end
    end
    return l
end

# EigenSpaceDecomposition

function eigen_decomposition!(M::Matrix{T}) where {T<:FiniteFields.GF}
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
        push!(eigspace_ptrs, cd + dim)
    end
    @assert eigspace_ptrs[end] == size(M, 1) + 1 "Matrix does not split over $T"
    return M, eigspace_ptrs
end

mutable struct EigenSpaceDecomposition{T<:FiniteFields.GF}
    basis::Matrix{T}
    eigspace_ptrs::Vector{Int}

    function EigenSpaceDecomposition(
        basis::Matrix{T},
        eigspace_ptrs::AbstractVector{<:Integer},
    ) where {T<:FiniteFields.GF}

        @assert eigspace_ptrs[1] == 1
        @assert eigspace_ptrs[end] == size(basis, 1) + 1
        return new{T}(basis, eigspace_ptrs)
    end
end

EigenSpaceDecomposition(M::Matrix{T}) where {T<:FiniteFields.GF} =
    EigenSpaceDecomposition(eigen_decomposition!(deepcopy(M))...)

function Base.show(
    io::IO,
    ::MIME"text/plain",
    esd::EigenSpaceDecomposition{T},
) where {T}
    println(io, tuple(diff(esd.eigspace_ptrs)...), "-splitting over ", T)
    print(io, esd.basis)
end

function Base.show(io::IO, esd::EigenSpaceDecomposition{T}) where {T}
    print(io, tuple(diff(esd.eigspace_ptrs)...), "-splitting over ", T)
end

Base.length(esd::EigenSpaceDecomposition) = length(esd.eigspace_ptrs) - 1

function Base.getindex(esd::EigenSpaceDecomposition, i::Int)
    @boundscheck 1 <= i <= length(esd)
    return esd.basis[esd.eigspace_ptrs[i]:esd.eigspace_ptrs[i+1]-1, :]
end

function Base.iterate(esd::EigenSpaceDecomposition, s = 1)
    s > length(esd) && return nothing
    first_last = esd.eigspace_ptrs[s]:esd.eigspace_ptrs[s+1]-1
    return (esd.basis[first_last, :], s + 1)
end

Base.eltype(esd::EigenSpaceDecomposition{T}) where {T} = Matrix{T}

LinearAlgebra.isdiag(esd::EigenSpaceDecomposition) =
    esd.eigspace_ptrs == 1:length(esd)+1

function refine(esd::EigenSpaceDecomposition{T}, M::Matrix{T}) where {T}
    nbasis = Array{T}(undef, 0, size(first(esd), 2))
    nptrs = [1]
    for (i, e) in enumerate(esd)
        if size(e, 1) > 1
            esd2, ptrs = eigen_decomposition!(e * @view M[:, _find_l(e)])
            nbasis = vcat(nbasis, esd2 * e)
            append!(nptrs, ptrs .+ (pop!(nptrs) - 1))
        else
            nbasis = vcat(nbasis, e)
            push!(nptrs, nptrs[end] + 1)
        end
    end
    return EigenSpaceDecomposition(nbasis, nptrs)
end
