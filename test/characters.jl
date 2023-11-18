@testset "Frobenius-Schur & projections" begin
    G = SmallPermGroups[8][4]
    irr = Characters.irreducible_characters(G)

    # quaternionic character
    χ = irr[1]
    @test collect(values(χ)) == [2, 0, -2, 0, 0]

    ι = Characters.frobenius_schur

    @test ι(χ) == -1

    @test deepcopy(χ) == χ
    @test deepcopy(χ) !== χ

    ψ = deepcopy(χ)

    @test ψ == χ
    @test ψ !== χ
    @test ψ == deepcopy(χ)

    @test Characters.conjugacy_classes(ψ) === Characters.conjugacy_classes(χ)
    @test Characters.table(ψ) === Characters.table(χ)

    @test collect(values(ψ)) == collect(values(χ))
    @test Characters.constituents(ψ) == Characters.constituents(χ)
    @test Characters.constituents(ψ) !== Characters.constituents(χ)

    @test hash(ψ) == hash(χ)

    @test 2ψ == 2χ
    @test 2ψ / 2 == χ
    @test 2ψ != χ

    @test PermutationGroups.degree(χ) == 2

    # @test size(SymbolicWedderburn.image_basis(χ), 1) == 4

    @test_throws AssertionError ι(2χ)
    a = 2χ
    @test sum(
        length(cc) * a(first(cc)^2) for cc in Characters.conjugacy_classes(a)
    ) // order(G) == 2ι(χ)

    # zero character
    @test zero(χ) == Characters.Character(Characters.table(χ), zeros(5))
    @test dot(χ, zero(χ)) == 0
    @test dot(zero(χ), zero(χ)) == 0
end

@testset "Characters io" begin
    G = PermGroup([perm"(2,3)(4,5)"])
    chars = Characters.irreducible_characters(G)

    @test sprint(show, chars[1]) == "χ₁"
    @test sprint(show, chars[2]) == "χ₂"
    @test sprint(show, chars[1] + 2chars[2]) == "χ₁ +2·χ₂"
    @test sprint(show, 2chars[1] - chars[2]) == "2·χ₁ -1·χ₂"
    @test sprint(show, -2chars[1] + 3chars[2]) == "-2·χ₁ +3·χ₂"

    @test sprint(show, MIME"text/plain"(), chars[1]) ==
          "Character over Cyclotomic{Rational{$Int}, SparseVector{Rational{$Int}, $Int}}\nχ₁"
end
