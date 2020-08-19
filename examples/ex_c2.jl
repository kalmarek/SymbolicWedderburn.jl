using SymbolicWedderburn
using PermutationGroups
using Cyclotomics


C2 = PermGroup([perm"(1,2)"])
ccG = conjugacy_classes(C2)
chars = SymbolicWedderburn.characters_dixon(ccG)

Cyclotomics.Cyclotomic{T, V}(a::R) where {T, V, R<:Real} = Cyclotomics.Cyclotomic{T, V}(1, [a])

Cyclotomics.Cyclotomic{T}(α::Cyclotomics.Cyclotomic) where {T} =
           Cyclotomics.Cyclotomic(conductor(α), convert.(T, α.coeffs))


function central_projection(chi::SymbolicWedderburn.AbstractClassFunction{T}) where T
    ccls = conjugacy_classes(chi)
    d = degree(one(first(first(ccls))))

    result = zeros(T, d, d)

    for cc in ccls
        val = chi(first(cc))
        for g in cc
            for i in 1:d
                result[i, i^g] += val
            end
        end
    end
    deg = degree(chi)
    ordG = sum(length, ccls)

    return deg//ordG, result
end

c1, U1 = central_projection(chars[1])
c2, U2 = central_projection(chars[2])

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




# @constraint m  sum(x.^2) - t in SOSCone() symmetry_group = C2
using MultivariateBases

mvec = monomials(x, 0:1)
dict = Dict(m => i for (i, m) in enumerate(mvec))

using MultivariatePolynomials
const MP = MultivariatePolynomials

function Base.:^(m::MP.AbstractMonomial, p::PermutationGroups.Perm)
    return prod(MP.variables(m).^(MP.exponents(m)^p))
end

G = PermGroup([Perm([dict[m] for m in mvec.^Ref(g)]) for g in gens(C2)])
ccG = conjugacy_classes(G)
chars = SymbolicWedderburn.characters_dixon(ccG)

c1, U1 = central_projection(chars[1])
c2, U2 = central_projection(chars[2])

#Int64 will not work for all Groups
R1, ids1 = SymbolicWedderburn.row_echelon_form(Int64.(U1))
R2, ids2 = SymbolicWedderburn.row_echelon_form(Int64.(U2))

msym = SOSModel(Mosek.Optimizer)
@variable   msym t
@objective  msym Max t
@variable   msym sos1 SOSPoly(FixedPolynomialBasis(R1[1:length(ids1),:]*mvec))
@variable   msym sos2 SOSPoly(FixedPolynomialBasis(R2[1:length(ids2),:]*mvec))
@constraint msym sum(x.^2) - t == sos1 + sos2
optimize!(msym)

# do we scale or not? 

# if we have a complex character, we need to find its conjugate
# 
