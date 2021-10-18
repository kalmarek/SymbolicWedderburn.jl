function zerotol!(M::AbstractArray; atol=eps(eltype(M)))
    @inbounds for i in eachindex(M)
        if abs(M[i]) < atol
            M[i] = zero(M[i])
        end
    end
    return M
end
