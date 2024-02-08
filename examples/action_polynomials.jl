using DynamicPolynomials
const DP = DynamicPolynomials
using GroupsCore

const SW = SymbolicWedderburn

# Defining action on polynomials by acting on terms and monomials:
function SW.action(a::SW.Action, el::GroupElement, poly::DP.AbstractPolynomial)
    return sum(SW.action(a, el, term) for term in DP.terms(poly))
end

function SW.action(a::SW.Action, el::GroupElement, term::DP.AbstractTerm)
    return DP.coefficient(term) * SW.action(a, el, DP.monomial(term))
end

struct VariablePermutation{V} <: SW.ByPermutations
    variables::V
end

function SW.action(
    a::VariablePermutation,
    g::AP.AbstractPermutation,
    m::Monomial,
)
    v = a.variables
    return m(v => SW.action(a, g, v))
end

# this is a general linear action that can be induced
# from the action on monomials
abstract type OnMonomials <: SW.ByLinearTransformation end

function SW.decompose(
    k::DP.AbstractPolynomialLike,
    hom::SW.InducedActionHomomorphism,
)
    # correct only if basis(hom) == monomials
    I = SW._int_type(hom)
    indcs = I[hom[mono] for mono in DP.monomials(k)]
    coeffs = DP.coefficients(k)

    return indcs, coeffs
end
