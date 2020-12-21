using SymbolicWedderburn
using PermutationGroups
using Cyclotomics

using DynamicPolynomials
using SumOfSquares
using SCS

@polyvar x y z

motzkin = x^4*y^2 + y^4*x^2 - 3*x^2*y^2 + 1

m = let m = SOSModel(optimizer_with_attributes(SCS.Optimizer, "eps"=>1e-5, "acceleration_lookback"=>10))

    @variable m t
    @objective m Max t
    @variable m sos SOSPoly(monomials([x, y], 0:5))
    @constraint m motzkin - t == sos
    optimize!(m)
    @info termination_status(m)
    m
end

include("action_polynomials.jl")

G = PermGroup([perm"(1,2)"])
basis = monomials([x,y], 0:3)

R = SymbolicWedderburn.symmetry_adapted_basis(G, basis)

msym = let R=R, msym = SOSModel(optimizer_with_attributes(SCS.Optimizer,
    "eps"=>3e-11,
    "max_iters"=>100_000,
    "acceleration_lookback"=>0,
    "alpha"=>1.95,
    ))

    @variable   msym t
    @objective  msym Max t

    @variable   msym sos1 SOSPoly(FixedPolynomialBasis(R[1]*basis))
    @variable   msym sos2 SOSPoly(FixedPolynomialBasis(R[2]*basis))

    @variable   msym sos3 SOSPoly(monomials([x,y], 0:1))

    @constraint msym motzkin - t + 0.125 == sos1 + sos2 + (2 - x^2 - y^2)*sos3
    @constraint msym t <= -0.05
    optimize!(msym)
    @info termination_status(msym)
    msym
end
