using Pkg; Pkg.activate(@__DIR__)
using Revise

using SymbolicWedderburn
using DynamicPolynomials
const MP = DynamicPolynomials.MP

using GroupsCore
include(joinpath(pkgdir(GroupsCore), "test", "cyclic.jl"))

G = CyclicGroup(2)
x,y = @polyvar x y
# (x, y) → (x+y, x-y)

struct MyAction <: SymbolicWedderburn.ByLinearTransformation end
# x -> linear combination of basis elements;

SymbolicWedderburn.coeff_type(m::MyAction) = Float64
function SymbolicWedderburn.action(a::MyAction, g::CyclicGroupElement, m::Monomial)
    isone(g) && return m
    x, y = variables(m)
    return m(x => (x+y)/sqrt(2), y => (x-y)/sqrt(2))
end

function SymbolicWedderburn.decompose(k::AbstractPolynomial, hom::SymbolicWedderburn.InducedActionHomomorphism)
    # correct only if features(hom) == monomials

    indcs = [hom[m] for m in monomials(k)]
    coeffs = coefficients(k)

    return indcs, coeffs
end



using Test

let g = gens(G, 1)

    basis = monomials([x,y], 0:4)
    k = SymbolicWedderburn.action(MyAction(), g, basis[2])
    ehom = SymbolicWedderburn.ExtensionHomomorphism(MyAction(), basis)
    idcs, vals = SymbolicWedderburn.decompose(k, ehom)
    @test sum(c*basis[i] for (c,i) in zip(vals, idcs)) == k

    for b in basis
        k = SymbolicWedderburn.action(MyAction(), g, b)
        idcs, vals = SymbolicWedderburn.decompose(k, ehom)
        @test sum(c*basis[i] for (c,i) in zip(vals, idcs)) == k
    end
end

# let us check that we get an actual matrix representation of G:

@testset "induced Matrix Representation" begin
    g = gens(G, 1)
    basis = monomials([x,y], 0:4)
    ehom = SymbolicWedderburn.ExtensionHomomorphism(MyAction(), basis)
    m = droptol!(SymbolicWedderburn.induce(MyAction(), ehom, g), 1e-30)
    @test eltype(m) == Float64

    @test m^2 ≈ one(m)

    B = SymbolicWedderburn.symmetry_adapted_basis(G, basis, MyAction());
    @test sum(first∘size, B) == length(basis)
end
