include(joinpath(pkgdir(GroupsCore), "test", "cyclic.jl"))

G = CyclicGroup(2)
x, y = @polyvar x y
# (x, y) → (x+y, x-y)

struct By90Rotation <: SW.ByLinearTransformation end
# x -> linear combination of basis elements;

SW.coeff_type(::By90Rotation) = Float64
function SW.action(::By90Rotation, g::CyclicGroupElement, m::AbstractMonomial)
    isone(g) && return m
    x, y = variables(m)
    return m(x => (x + y) / sqrt(2), y => (x - y) / sqrt(2)) # we need the action to be orthogonal
end

function SW.decompose(
    k::DP.AbstractPolynomialLike,
    hom::SW.InducedActionHomomorphism,
)
    # correct only if basis(hom) == monomials
    I = SW._int_type(hom)
    indcs = I[hom[mono] for mono in DP.monomials(k)]
    coeffs = DP.coefficients(k)

    return indcs, coeffs
end

using Test

@testset "Decompose in basis" begin
    g = gens(G, 1)

    basis = SA.FixedBasis(monomials([x, y], 0:4), SA.DiracMStructure(*))
    k = SW.action(By90Rotation(), g, basis[2])
    ehom = SW.ExtensionHomomorphism(By90Rotation(), basis)
    idcs, vals = SW.decompose(k, ehom)
    @test sum(c * basis[i] for (c, i) in zip(vals, idcs)) == k

    for b in basis
        k = SW.action(By90Rotation(), g, b)
        idcs, vals = SW.decompose(k, ehom)
        @test sum(c * basis[i] for (c, i) in zip(vals, idcs)) == k
    end
end

# let us check that we get an actual matrix representation of G:

@testset "induced Matrix Representation" begin
    g = gens(G, 1)
    monomial_basis =
        SA.FixedBasis(monomials([x, y], 0:4), SA.DiracMStructure(*))
    ehom = SW.ExtensionHomomorphism(By90Rotation(), monomial_basis)
    m = droptol!(SW.induce(By90Rotation(), ehom, g), 1e-15)
    @test eltype(m) == Float64

    @test !(m ≈ one(m))
    @test m^2 ≈ one(m)

    ssimple_basis =
        SW.symmetry_adapted_basis(Float64, G, By90Rotation(), monomial_basis)
    degs = SW.degree.(ssimple_basis)
    mlts = multiplicity.(ssimple_basis)
    @test sum(d * m for (d, m) in zip(degs, mlts)) == length(monomial_basis)
    @test sum(first ∘ size, ssimple_basis) == length(monomial_basis)
end
