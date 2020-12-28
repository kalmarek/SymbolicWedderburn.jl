using DynamicPolynomials
using MultivariatePolynomials

function SymbolicWedderburn.InducingHomomorphism(basis::AbstractVector{<:AbstractMonomial})
    basis_exps = Vector{Vector{Int}}(undef, length(basis))
    basis_dict = Dict{Vector{Int},Int}()
    sizehint!(basis_dict, length(basis))

    for (i, b) in enumerate(basis)
        e = MultivariatePolynomials.exponents(b) # so that we allocate exponents only once
        basis_exps[i] = e
        basis_dict[e] = i
    end

    return SymbolicWedderburn.InducingHomomorphism(basis_exps, basis_dict)
end
