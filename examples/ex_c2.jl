# Dense example
using DynamicPolynomials
using SumOfSquares
using MosekTools

@polyvar x[1:2]

m = SOSModel(Mosek.Optimizer)
@variable m t
@objective m Max t
@variable m sos SOSPoly([1, x[1], x[2]])
@constraint m sum(x.^2) - t == sos
optimize!(m)

# Exploit Symmetry
using SymbolicWedderburn
using PermutationGroups
using Cyclotomics

msym = SOSModel(Mosek.Optimizer)
@variable   msym t
@objective  msym Max t

C2 = PermGroup([perm"(1,2)"])

# @constraint msym sum(x.^2) - t in SOSCone() symmetry_group = C2

using MultivariateBases
using MultivariatePolynomials
const MP = MultivariatePolynomials

mvec = monomials(x, 0:1)
dict = Dict(m => i for (i, m) in enumerate(mvec))

function Base.:^(m::MP.AbstractMonomial, p::PermutationGroups.Perm)
    return prod(MP.variables(m).^(MP.exponents(m)^p))
end

G = PermGroup([Perm([dict[m] for m in mvec.^Ref(g)]) for g in gens(C2)])
ccG = conjugacy_classes(G)
chars = SymbolicWedderburn.characters_dixon(ccG)

c1, U1 = SymbolicWedderburn.central_projection(chars[1])
c2, U2 = SymbolicWedderburn.central_projection(chars[2])

#Int64 will not work for all Groups
R1, ids1 = SymbolicWedderburn.row_echelon_form(Int64.(U1))
R2, ids2 = SymbolicWedderburn.row_echelon_form(Int64.(U2))

@variable msym sos1 SOSPoly(FixedPolynomialBasis(R1[1:length(ids1),:]*mvec))
@variable msym sos2 SOSPoly(FixedPolynomialBasis(R2[1:length(ids2),:]*mvec))
@constraint msym sum(x.^2) - t == sos1 + sos2
optimize!(msym)

# do we scale or not? 

# if we have a complex character, we need to find its conjugate
# 
