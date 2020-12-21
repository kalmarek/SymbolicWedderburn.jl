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
    @constraint m sum(x) + sum(x .^ 2) - t == sos
    optimize!(m)
end

# @constraint m  sum(x.^2) - t in SOSCone() symmetry_group = C4

include("action_polynomials.jl")

# include(joinpath(@__DIR__, "..", "test", "smallgroups.jl"));
# G = SmallPermGroups[4][1]

G = PermGroup([perm"(1,2,3,4)"])
# G = PermGroup([perm"(1,2)", perm"(1,2,3,4)"])

basis = monomials(x, 0:1)

R = SymbolicWedderburn.symmetry_adapted_basis(G, basis)

msym = let msym = SOSModel(SCS.Optimizer), basis = basis, R = R
    @variable msym t
    @objective msym Max t

    @variable msym sos1 SOSPoly(FixedPolynomialBasis(R[1] * basis))
    @variable msym sos2 SOSPoly(FixedPolynomialBasis(R[2] * basis))
    @variable msym sos3 SOSPoly(FixedPolynomialBasis(R[3] * basis))

    @constraint msym sum(x) + sum(x .^ 2) - t == sos1 + sos2 + sos3
    optimize!(msym)
    @info termination_status(msym)
    msym
end

value.(msym[:sos1].Q)

value.(msym[:sos2].Q)

value.(msym[:sos3].Q)

