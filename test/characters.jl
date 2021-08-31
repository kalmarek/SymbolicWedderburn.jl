@testset "Frobenius-Schur & projections" begin
    G = SmallPermGroups[8][4]
    irr = SymbolicWedderburn.irreducible_characters(G)

    # quaternionic character
    χ = irr[1]
    @test collect(values(χ)) == [2, 0, 0, -2, 0]

    ι = SymbolicWedderburn.frobenius_schur

    @test ι(χ) == -1

    @test deepcopy(χ) == χ
    @test deepcopy(χ) !== χ

    ψ = deepcopy(χ)

    @test ψ == χ
    @test ψ !== χ
    @test ψ == deepcopy(χ)

    @test conjugacy_classes(ψ) === conjugacy_classes(χ)
    @test SymbolicWedderburn.table(ψ) === SymbolicWedderburn.table(χ)

    @test collect(values(ψ)) == collect(values(χ))
    @test SymbolicWedderburn.constituents(ψ) == SymbolicWedderburn.constituents(χ)
    @test SymbolicWedderburn.constituents(ψ) !== SymbolicWedderburn.constituents(χ)

    @test hash(ψ) == hash(χ)

    @test 2ψ == 2χ
    @test 2ψ/2 == χ
    @test 2ψ != χ

    @test PermutationGroups.degree(χ) == 2

    @test size(SymbolicWedderburn.image_basis(χ), 1) == 4

    @test_throws AssertionError ι(2χ)
    a = 2χ
    @test sum(length(cc)*a(first(cc)^2) for cc in conjugacy_classes(a))//order(G) == 2ι(χ)
end
