permuted(exponent::AbstractVector, p::Perm) =
    [exponent[i^p] for i in eachindex(exponent)]
