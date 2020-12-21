using SymbolicWedderburn
using PermutationGroups
using Cyclotomics

using DynamicPolynomials
using SumOfSquares
using SCS

@polyvar x[1:4]

let m = SOSModel(SCS.Optimizer)
    @variable m t
    @objective m Max t
    @variable m sos SOSPoly([1; x])
    @constraint m sum(x) + sum(x.^2) - t == sos
    optimize!(m)
end

# @constraint m  sum(x.^2) - t in SOSCone() symmetry_group = C2
using DynamicPolynomials
using MultivariatePolynomials
const MP = MultivariatePolynomials
using MultivariateBases


include("action_polynomials.jl")


include(joinpath(@__DIR__, "..", "test", "smallgroups.jl"));
G = SmallPermGroups[4][1]
# G = PermGroup([perm"(1,2)", perm"(1,2,3,4)"])
chars_vars = SymbolicWedderburn.characters_dixon(G)

mvec = monomials(x, 0:1)

chars_mvec = let chars = chars_vars, basis = mvec

    @assert all(χ.inv_of == first(chars).inv_of for χ in chars)

    induced_action = SymbolicWedderburn.OnPoints(basis)
    ccG_large = induced_action.(conjugacy_classes(first(chars)))

    # double check:
    let ccls = ccG_large, large_gens = induced_action.(gens(G))
        G_large = PermGroup(large_gens)
        ccG_large = conjugacy_classes(G_large)
        @assert all(Set.(collect.(ccG_large)) .== Set.(collect.(ccls)))
    end

    [SymbolicWedderburn.Character(values(χ), χ.inv_of, ccG_large) for χ in chars]
end

vr_chars = SymbolicWedderburn.real_vchars(chars_mvec)

U = filter!(x->all(!iszero, x), [SymbolicWedderburn.matrix_projection(χ) for χ in vr_chars])

R = map(U) do c_u
    u = last(c_u)
    if all(isreal, u)
        image_coeffs, pivots = SymbolicWedderburn.row_echelon_form(float.(u))
        dim = length(pivots)
        image_coeffs[1:dim, :]
    else
        throw("Not Implemented")
    end
end

R = filter!(!iszero ∘ first ∘ size, R)

msym = let msym = SOSModel(SCS.Optimizer)
    @variable   msym t
    @objective  msym Max t

    @variable   msym sos1 SOSPoly(FixedPolynomialBasis(R[1]*mvec))
    @variable   msym sos2 SOSPoly(FixedPolynomialBasis(R[2]*mvec))
    @variable   msym sos3 SOSPoly(FixedPolynomialBasis(R[3]*mvec))
    # @variable   msym sos4
    @constraint msym sum(x) + sum(x.^2) - t == sos1 + sos2 + sos3
    optimize!(msym)
    msym
end

termination_status(msym)

value.(msym[:sos1].Q)

value.(msym[:sos2].Q)

value.(msym[:sos3].Q)

