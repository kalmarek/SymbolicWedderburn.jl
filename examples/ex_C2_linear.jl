using SymbolicWedderburn
using DynamicPolynomials
using SparseArrays

using GroupsCore
include(joinpath(pkgdir(GroupsCore), "test", "cyclic.jl"))

G = CyclicGroup(2)
x,y = @polyvar x y
# (x, y) → (x+y, x-y)

struct By90Rotation <: SymbolicWedderburn.ByLinearTransformation end
# x -> linear combination of basis elements;

SymbolicWedderburn.coeff_type(::By90Rotation) = Float64
function SymbolicWedderburn.action(::By90Rotation, g::CyclicGroupElement, m::Monomial)
    isone(g) && return m
    x, y = variables(m)
    return m(x => (x+y)/sqrt(2), y => (x-y)/sqrt(2)) # we need the action to be orthogonal
end

function SymbolicWedderburn.decompose(k::AbstractPolynomial, hom::SymbolicWedderburn.InducedActionHomomorphism)
    # correct only if features(hom) == monomials

    indcs = [hom[m] for m in monomials(k)]
    coeffs = coefficients(k)

    return indcs, coeffs
end



using Test

@testset "Decompose in basis" begin
    g = gens(G, 1)

    basis = monomials([x,y], 0:4)
    k = SymbolicWedderburn.action(By90Rotation(), g, basis[2])
    ehom = SymbolicWedderburn.ExtensionHomomorphism(By90Rotation(), basis)
    idcs, vals = SymbolicWedderburn.decompose(k, ehom)
    @test sum(c*basis[i] for (c,i) in zip(vals, idcs)) == k

    for b in basis
        k = SymbolicWedderburn.action(By90Rotation(), g, b)
        idcs, vals = SymbolicWedderburn.decompose(k, ehom)
        @test sum(c*basis[i] for (c,i) in zip(vals, idcs)) == k
    end
end

# let us check that we get an actual matrix representation of G:

@testset "induced Matrix Representation" begin
    g = gens(G, 1)
    basis = monomials([x,y], 0:4)
    ehom = SymbolicWedderburn.ExtensionHomomorphism(By90Rotation(), basis)
    m = droptol!(SymbolicWedderburn.induce(By90Rotation(), ehom, g), 1e-15)
    @test eltype(m) == Float64

    @test !(m ≈ one(m))
    @test m^2 ≈ one(m)

    ssimple_basis = SymbolicWedderburn.symmetry_adapted_basis(Float64, G, By90Rotation(), basis);
    degs = SymbolicWedderburn.degree.(ssimple_basis)
    mlts = SymbolicWedderburn.multiplicity.(ssimple_basis)
    @test sum(d*m for (d,m) in zip(degs, mlts)) == length(basis)
    @test sum(first∘size∘SymbolicWedderburn.basis, ssimple_basis) == length(basis)
end
