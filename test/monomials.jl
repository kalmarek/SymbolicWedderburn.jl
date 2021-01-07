#= 
function Base.:^(m::MultivariatePolynomials.AbstractMonomialLike, g::Perm{Int64})

# This function does not work at the moment, since DynamicPolynomials sometimes reduces the number of variables (see https://github.com/JuliaAlgebra/DynamicPolynomials.jl/issues/79 )

return prod(variables(m).^(exponents(m)^g))
end

function Base.:^(t::MultivariatePolynomials.AbstractTermLike, g::Perm{Int64})
return coefficient(t)*monomial(t)^g
end

function Base.:^(p::MultivariatePolynomials.AbstractPolynomialLike, g::Perm{Int64})
return sum([t^g for t in terms(p)])    
end

function test_invariant(G, poly)
@show poly
@show variables(poly)
# should probably only iterate over a subgroup
for g in Iterators.drop(G, 0)
@show(g)
@test poly^g == poly
end
end
=#

@testset "Symmetry on PolyVar" begin
    @polyvar x[1:2]
    @polyvar y[1:2]
    basis = [x; y]
    G = PermGroup([perm"(1,2,3,4)"])
    R = symmetry_adapted_basis_float(G)
    new_bases = map(R) do Ri
        FixedPolynomialBasis(Ri * basis)
    end

    @test length.(new_bases) == [1, 1, 2]
    @test new_bases[1].polynomials[1] == x[1] - x[2] + y[1] - y[2]
    @test new_bases[2].polynomials[1] == x[1] + x[2] + y[1] + y[2]
    @test new_bases[3].polynomials[1] == x[1] - y[1]
    @test new_bases[3].polynomials[2] == x[2] - y[2]

end


function SymbolicWedderburn.InducingHomomorphism(basis::AbstractVector{<:AbstractMonomial})
    basis_exps = Vector{Vector{Int}}(undef, length(basis))
    basis_dict = Dict{Vector{Int},Int}()
    sizehint!(basis_dict, length(basis))

    for (i, b) in enumerate(basis)
        e = MultivariatePolynomials.exponents(b) # so that we allocate exponents only once
        basis_exps[i] = e
        basis_dict[e] = i
    end

    return SymbolicWedderburn.InducingHomomorphism(basis_exps, basis_dict)
end

@testset "Symmetry on Monomials" begin
    @polyvar x[1:2]
    basis = monomials(x, 0:2)
    G = PermGroup([perm"(1,2)"])
    R = symmetry_adapted_basis_float(G, basis)
    new_bases = map(R) do Ri
        FixedPolynomialBasis(Ri * basis)
    end

    @test length.(new_bases) == [2, 4]
    @test new_bases[1].polynomials[1] == x[1]^2 - x[2]^2
    @test new_bases[1].polynomials[2]== x[1] - x[2]
    @test new_bases[2].polynomials[1] ==  x[1]^2 + x[2]^2
    @test new_bases[2].polynomials[2] ==  x[1] * x[2]
    @test new_bases[2].polynomials[3] ==  x[1] + x[2]
    @test new_bases[2].polynomials[4] == 1
end