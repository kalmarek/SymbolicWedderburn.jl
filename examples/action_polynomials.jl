using DynamicPolynomials
const DPoly = DynamicPolynomials
using GroupsCore
using PermutationGroups

struct VariablePermutation <: SymbolicWedderburn.ByPermutations end

function SymbolicWedderburn.action(
    a::VariablePermutation,
    g::PermutationGroups.AbstractPerm,
    m::Monomial,
)
    v = variables(m)
    return m(v => SymbolicWedderburn.action(a, g, v))
end

function SymbolicWedderburn.action(
    ::VariablePermutation,
    g::PermutationGroups.AbstractPerm,
    v::AbstractVector,
)
    return map(i -> v[i^g], eachindex(v))
end

abstract type OnMonomials <: SymbolicWedderburn.ByLinearTransformation end

function SymbolicWedderburn.action(
    a::Union{VariablePermutation,OnMonomials},
    el::GroupElement,
    term::DPoly.AbstractTerm,
)
    return DPoly.coefficient(term) *
           SymbolicWedderburn.action(a, el, DPoly.monomial(term))
end
function SymbolicWedderburn.action(
    a::Union{VariablePermutation,OnMonomials},
    el::GroupElement,
    poly::DPoly.AbstractPolynomial,
)
    return sum([
        SymbolicWedderburn.action(a, el, term) for term in DPoly.terms(poly)
    ])
end

function SymbolicWedderburn.decompose(
    k::DPoly.AbstractPolynomialLike,
    hom::SymbolicWedderburn.InducedActionHomomorphism,
)
    # correct only if basis(hom) == monomials
    I = SymbolicWedderburn._int_type(hom)
    indcs = I[hom[mono] for mono in DPoly.monomials(k)]
    coeffs = DPoly.coefficients(k)

    return indcs, coeffs
end
