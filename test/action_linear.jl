using DynamicPolynomials

using GroupsCore
include(joinpath(pathof(GroupsCore), "..", "..", "test", "cyclic.jl"))

x,y = @polyvar x y
struct MyAction <: SymbolicWedderburn.ByLinearTransformation end
# (x, y) → (x+y, x-y)

SymbolicWedderburn.coeff_type(m::MyAction) = Float64
function SymbolicWedderburn.action(a::MyAction, g::CyclicGroupElement, m::Monomial)
    isone(g) && return m
    x, y = variables(m)
    # we need the action to be orthogonal
    return m(x => (x+y)/sqrt(2), y => (x-y)/sqrt(2))
end

function SymbolicWedderburn.decompose(k::AbstractPolynomial, hom::SymbolicWedderburn.InducedActionHomomorphism)
    # correct only if features(hom) == monomials

    indcs = [hom[m] for m in monomials(k)]
    coeffs = coefficients(k)

    return indcs, coeffs
end

@testset "Linear Actions" begin

    G = CyclicGroup(2)
    basis = monomials([x,y], 0:4)
    ehom = SymbolicWedderburn.ExtensionHomomorphism(MyAction(), basis)

    @testset "decomposition in basis" begin
        g = gens(G, 1)

        k = SymbolicWedderburn.action(MyAction(), g, basis[2])
        idcs, vals = SymbolicWedderburn.decompose(k, ehom)
        @test sum(c*basis[i] for (c,i) in zip(vals, idcs)) == k

        for b in basis
            k = SymbolicWedderburn.action(MyAction(), g, b)
            idcs, vals = SymbolicWedderburn.decompose(k, ehom)
            @test sum(c*basis[i] for (c,i) in zip(vals, idcs)) == k
        end
    end

    @testset "induced Matrix Representation" begin
        g = gens(G, 1)
        M = SymbolicWedderburn.induce(MyAction(), ehom, g)

        @test M isa AbstractMatrix{Float64}
        @test size(M) == (length(basis), length(basis))

        @test !(M ≈ one(M))
        @test M^2 ≈ one(M)

        B = SymbolicWedderburn.symmetry_adapted_basis(G, basis, MyAction());
        @test sum(first ∘ size ∘ SymbolicWedderburn.basis, B) == length(basis)
    end
end
