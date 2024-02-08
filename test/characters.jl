@testset "Frobenius-Schur & projections" begin
    G = SmallPermGroups[8][4]
    tbl = Characters.CharacterTable(Rational{Int}, G)

    # quaternionic character
    χ = Characters.Character{Rational{Int}}(tbl, 5)
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
    @test Characters.multiplicities(ψ) == Characters.multiplicities(χ)
    @test Characters.multiplicities(ψ) !== Characters.multiplicities(χ)

    @test hash(ψ) == hash(χ)

    @test 2ψ == 2χ
    @test 2ψ / 2 == χ
    @test 2ψ != χ

    @test AP.degree(χ) == 2

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

    @testset "characters multiplication: Sym(3)" begin
        G = PermGroup(perm"(1,2,3)", perm"(1,2)")
        tbl = Characters.CharacterTable(G)
        χ = Characters.irreducible_characters(tbl)

        # χ[1] - the trivial character
        # χ[2] - the alternating character
        # χ[3] - the non-trivial, degree-2 character

        @test all(χ[1] * ψ == ψ for ψ in χ)
        @test χ[2]^2 == χ[1]
        @test χ[2] * χ[3] == χ[3] * χ[2] == χ[3]

        @test χ[3]^2 == χ[1] + χ[2] + χ[3]
    end

    @testset "characters multiplication: Sym(4)" begin
        G = PermGroup(perm"(1,2,3,4)", perm"(1,2)")
        tbl = Characters.CharacterTable(G)
        χ = Characters.irreducible_characters(tbl)

        # χ[1] - the trivial character
        # χ[2] - the alternating character
        # χ[3] - the non-trivial, degree-2 character
        # χ[4] - the non-trivial, degree-3 character
        # χ[5] - χ[4]*χ[2]

        @test all(χ[1] * ψ == ψ for ψ in χ)
        @test χ[2]^2 == χ[1]
        @test χ[2] * χ[3] == χ[3]
        @test χ[2] * χ[4] == χ[5]
        @test χ[2] * χ[5] == χ[4]

        @test χ[3]^2 == χ[1] + χ[2] + χ[3]
        @test χ[3] * χ[4] == χ[4] + χ[5]
        @test χ[3] * χ[5] == χ[4] + χ[5]

        @test χ[4]^2 == χ[1] + χ[3] + χ[4] + χ[5]
        @test χ[4] * χ[5] == χ[2] + χ[3] + χ[4] + χ[5]

        @test χ[5]^2 == χ[1] + χ[3] + χ[4] + χ[5]
    end
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
          "Character over cyclotomics (Rational{$Int})\nχ₁"

    χ = Characters.Character{Rational{BigInt}}(chars[1])
    @test sprint(show, MIME"text/plain"(), χ) ==
          "Character over rationals ($BigInt)\nχ₁"

    @test sprint(show, MIME"text/plain"(), Characters.table(χ)) isa
          AbstractString
end
