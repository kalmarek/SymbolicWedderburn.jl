using DynamicPolynomials
const  MPoly = DynamicPolynomials.MultivariatePolynomials
using GroupsCore
using PermutationGroups

struct VariablePermutation <: SymbolicWedderburn.ByPermutations end

function SymbolicWedderburn.action(a::VariablePermutation, g::PermutationGroups.AbstractPerm, m::Monomial)
    v = variables(m)
    return m(v => SymbolicWedderburn.action(a, g, v))
end

function SymbolicWedderburn.action(::VariablePermutation, g::PermutationGroups.AbstractPerm, v::AbstractVector)
    return map(i -> v[i^g], eachindex(v))
end

abstract type OnMonomials <: SymbolicWedderburn.ByLinearTransformation end

function SymbolicWedderburn.action(
    a::Union{VariablePermutation,OnMonomials},
    el::GroupElement,
    term::MPoly.AbstractTerm,
)
    return MPoly.coefficient(term) * SymbolicWedderburn.action(a, el, MPoly.monomial(term))
end
function SymbolicWedderburn.action(
    a::Union{VariablePermutation,OnMonomials},
    el::GroupElement,
    poly::MPoly.AbstractPolynomial,
)
    return sum([SymbolicWedderburn.action(a, el, term) for term in MPoly.terms(poly)])
end

function SymbolicWedderburn.decompose(
    k::MPoly.AbstractPolynomialLike,
    hom::SymbolicWedderburn.InducedActionHomomorphism,
)
    # correct only if basis(hom) == monomials
    I = SymbolicWedderburn._int_type(hom)
    indcs = I[hom[mono] for mono in MPoly.monomials(k)]
    coeffs = MPoly.coefficients(k)

    return indcs, coeffs
end
