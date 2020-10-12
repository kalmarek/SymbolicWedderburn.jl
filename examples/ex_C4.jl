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

function SymbolicWedderburn.OnPoints(basis::Union{MonomialVector, AbstractVector{<:Monomial}})
    basis_exps = Vector{Vector{Int}}(undef, length(basis))
    basis_dict = Dict{Vector{Int}, Int}()
    sizehint!(basis_dict, length(basis))

    for (i, b) in enumerate(basis)
        e = MP.exponents(b) # so that we allocate exponents only once
        basis_exps[i] = e
        basis_dict[e] = i
    end

    return SymbolicWedderburn.OnPoints(basis_exps, basis_dict)
end

include(joinpath(@__DIR__, "..", "test", "smallgroups.jl"));
G = SmallPermGroups[4][1]
chars_vars = SymbolicWedderburn.characters_dixon(G)

mvec = reverse(monomials(x, 0:1))

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

vr_chars = SymbolicWedderburn.real_vchars(chars_large)
U = [SymbolicWedderburn.matrix_projection(χ) for χ in vr_chars]

R = map(U) do c_u
    u = last(c_u)
    if all(isreal, u)
        SymbolicWedderburn.row_echelon_form(float.(u))
    else
        throw("Not Implemented")
    end
end

msym = let msym = SOSModel(SCS.Optimizer)
    @variable   msym t
    @objective  msym Max t
    @variable   msym sos1 SOSPoly(FixedPolynomialBasis(first(R[1])[1:length(last(R[1])),:]*mvec))
    @variable   msym sos2 SOSPoly(FixedPolynomialBasis(first(R[2])[1:length(last(R[2])),:]*mvec))
    @variable   msym sos3 SOSPoly(FixedPolynomialBasis(first(R[3])[1:length(last(R[3])),:]*mvec))
    @variable   msym sos4 SOSPoly(FixedPolynomialBasis(first(R[4])[1:length(last(R[4])),:]*mvec))
    @constraint msym sum(x) + sum(x.^2) - t == sos1 + sos2 + sos3 + sos4
    optimize!(msym)
    msym
end
