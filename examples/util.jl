using SparseArrays

function invariant_constraint!(
    M_orb::AbstractMatrix,
    M::Matrix{<:Integer},
    invariant_vec::SparseVector,
)
    M_orb .= zero(eltype(M_orb))
    for i in eachindex(M)
        if M[i] ∈ SparseArrays.nonzeroinds(invariant_vec)
            M_orb[i] += invariant_vec[M[i]]
        end
    end
    return M_orb
end

function invariant_constraint!(
    M_orb::AbstractMatrix,
    M::Matrix,
    invariant_vec::AbstractVector,
)
    M_orb .= zero(eltype(M_orb))
    for i in eachindex(M)
        M_orb[i] += invariant_vec[M[i]]
    end
    return M_orb
end
