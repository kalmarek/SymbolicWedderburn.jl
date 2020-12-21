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
    m
end

include("action_polynomials.jl")


include(joinpath(@__DIR__, "..", "test", "smallgroups.jl"));
G = PermGroup([perm"(1,2)"])
chars_vars = SymbolicWedderburn.characters_dixon(G)

mvec = monomials([x,y], 0:3)

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

msym = let msym = SOSModel(optimizer_with_attributes(SCS.Optimizer,
    "eps"=>1e-12,
    "max_iters"=>100_000,
    "acceleration_lookback"=>0,
    "alpha"=>1.95,
    ))

    @variable   msym t
    @objective  msym Max t

    @variable   msym sos1 SOSPoly(FixedPolynomialBasis(R[1]*mvec))
    @variable   msym sos2 SOSPoly(FixedPolynomialBasis(R[2]*mvec))

    @variable   msym sos3 SOSPoly(monomials([x,y], 0:1))

    @constraint msym motzkin - t + 0.125 == sos1 + sos2 + (2 - x^2 - y^2)*sos3
    @constraint msym t <= -0.05
    optimize!(msym)
    msym
end


termination_status(msym)

