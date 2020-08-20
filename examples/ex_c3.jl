using SymbolicWedderburn
using PermutationGroups
using Cyclotomics


using DynamicPolynomials
using SumOfSquares
using SCS


C3 = PermGroup([perm"(1,2,3)"])
@polyvar x[1:3]

# @constraint m  sum(x.^2) - t in SOSCone() symmetry_group = C2
using MultivariateBases

mvec = monomials(x, 0:1)
dict = Dict(m => i for (i, m) in enumerate(mvec))

using MultivariatePolynomials
const MP = MultivariatePolynomials

function Base.:^(m::MP.AbstractMonomial, p::PermutationGroups.Perm)
    return prod(MP.variables(m).^(MP.exponents(m)^p))
end

G = PermGroup([Perm([dict[m] for m in mvec.^Ref(g)]) for g in gens(C3)])
ccG = conjugacy_classes(G)
chars = SymbolicWedderburn.characters_dixon(ccG)

c1, U1 = central_projection(chars[1])
c2, U2 = central_projection(chars[2])
c3, U3 = central_projection(chars[3])

@assert conj(U2) ==  U3

R1, ids1 = SymbolicWedderburn.row_echelon_form(float(U1))
R2, ids2 = SymbolicWedderburn.row_echelon_form(float(U2+U3))
R3, ids3 = SymbolicWedderburn.row_echelon_form(float(-E(4)*(U2-U3)))

msym = SOSModel(SCS.Optimizer)
@variable   msym t
@objective  msym Max t
@variable   msym sos1 SOSPoly(FixedPolynomialBasis(R1[1:length(ids1),:]*mvec))
@variable   msym sos2 SOSPoly(FixedPolynomialBasis(R2[1:length(ids2),:]*mvec))
@variable   msym sos3 SOSPoly(FixedPolynomialBasis(R3[1:length(ids3),:]*mvec))

@constraint msym sum(x.^2) - t == sos1 + sos2 + sos3
optimize!(msym)
#=  
    Variables n = 10, constraints m = 19
    Cones:	primal zero / dual free vars: 10
	        soc vars: 9, soc blks: 3  
=#

m = let 
    m = SOSModel(SCS.Optimizer)
    @variable m t
    @objective m Max t
    @constraint m  sum(x.^2) - t in SOSCone()
    optimize!(m)

    #=
    Variables n = 11, constraints m = 20
    Cones:	primal zero / dual free vars: 10
	        sd vars: 10, sd blks: 1
    =#
end


