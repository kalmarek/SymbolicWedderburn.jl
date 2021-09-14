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

function _findmax(f, A::AbstractArray)
    mval, midx = f(first(A)), firstindex(A)
    for idx in eachindex(A)
        if (v = f(A[idx])) > mval
            mval, midx = v, idx
        end
    end
    return mval, midx
end

function row_echelon_form!(A::AbstractMatrix{T}, atol=eps(T)*max(size(A)...)) where {T<:AbstractFloat}
    pivots = Int[]
    pos = 0
    sgn = 1.0
    for i = 1:size(A, 2)
        # find the largest entry in the column below the last pivot
        mval, midx = _findmax(abs, @view(A[pos+1:end, i]))
        # j = findnext(x -> !iszero(x), , pos+1)
        if mval < atol # the largest entry is below threshold so we zero everything in the column
            A[pos+1:end, i] .= zero(T)
            continue
        end
        j = midx+pos
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
        sgn *= sign(w)
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
    if sgn < 0
        A .*= sgn
    end
    return A, pivots
end

row_echelon_form(A::AbstractMatrix) = row_echelon_form!(deepcopy(A))
