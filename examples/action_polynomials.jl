using DynamicPolynomials
MP = DynamicPolynomials.MultivariatePolynomials
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
    term::MP.AbstractTerm,
)
    return MP.coefficient(term) * SymbolicWedderburn.action(a, el, MP.monomial(term))
end
function SymbolicWedderburn.action(
    a::Union{VariablePermutation,OnMonomials},
    el::GroupElement,
    poly::MP.AbstractPolynomial,
)
    return sum([SymbolicWedderburn.action(a, el, term) for term in MP.terms(poly)])
end

function SymbolicWedderburn.decompose(
    k::MP.AbstractPolynomialLike,
    hom::SymbolicWedderburn.InducedActionHomomorphism,
)
    # correct only if basis(hom) == monomials

    indcs = [hom[mono] for mono in MP.monomials(k)]
    coeffs = MP.coefficients(k)

    return indcs, coeffs
end
