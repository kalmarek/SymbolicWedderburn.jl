using DynamicPolynomials

using GroupsCore
include(joinpath(dirname(pathof(GroupsCore)), "..", "test", "cyclic.jl"))

x, y = @polyvar x y
struct By90Rotation <: SymbolicWedderburn.ByLinearTransformation end
# It would be clearer to make this action <: OnMonomials
# (x, y) → (x+y, x-y)/sqrt(2)

SymbolicWedderburn.coeff_type(::By90Rotation) = Float64
function SymbolicWedderburn.action(
    ::By90Rotation,
    g::CyclicGroupElement,
    m::AbstractMonomial,
)
    isone(g) && return m
    x, y = variables(m)
    return m(x => (x + y) / sqrt(2), y => (x - y) / sqrt(2)) # we need the action to be orthogonal
end

function SymbolicWedderburn.decompose(
    k::DynamicPolynomials.AbstractPolynomialLike,
    hom::SymbolicWedderburn.InducedActionHomomorphism,
)
    # correct only if basis(hom) == monomials
    I = SymbolicWedderburn._int_type(hom)
    indcs = I[hom[mono] for mono in DynamicPolynomials.monomials(k)]
    coeffs = DynamicPolynomials.coefficients(k)

    return indcs, coeffs
end

@testset "Linear Actions" begin
    G = CyclicGroup(2)
    monomial_basis = let monoms = monomials([x, y], 0:4)
        SA.FixedBasis(monoms, SA.DiracMStructure(*))
    end
    ehom =
        SymbolicWedderburn.ExtensionHomomorphism(By90Rotation(), monomial_basis)

    @testset "decomposition in basis" begin
        g = gens(G, 1)

        k = SymbolicWedderburn.action(By90Rotation(), g, monomial_basis[2])
        idcs, vals = SymbolicWedderburn.decompose(k, ehom)
        @test sum(c * monomial_basis[i] for (c, i) in zip(vals, idcs)) == k

        for b in monomial_basis
            k = SymbolicWedderburn.action(By90Rotation(), g, b)
            idcs, vals = SymbolicWedderburn.decompose(k, ehom)
            @test sum(c * monomial_basis[i] for (c, i) in zip(vals, idcs)) == k
        end
    end

    @testset "induced Matrix Representation" begin
        g = gens(G, 1)
        M = SymbolicWedderburn.induce(By90Rotation(), ehom, g)

        @test M isa AbstractMatrix{Float64}
        @test size(M) == (length(monomial_basis), length(monomial_basis))

        @test !(M ≈ one(M))
        @test M^2 ≈ one(M)

        sa_basis =
            symmetry_adapted_basis(Float64, G, By90Rotation(), monomial_basis)
        @test sum(first ∘ size, sa_basis) == length(monomial_basis)
    end
end
