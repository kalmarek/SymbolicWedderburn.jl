function generictest_dixon_Fp(G, p = Characters.dixon_prime(G))
    F = Characters.FiniteFields.GF{p}
    ccG = Characters.conjugacy_classes(G)
    Ns = [Characters.CMMatrix(ccG, i) for i in 1:length(ccG)]
    @test isdiag(Characters.common_esd(Ns, F))
    esd = Characters.common_esd(Ns, F)

    let W = esd.basis
        @test inv(W) isa Matrix{F}
        @test isone(inv(W) * W)
        iW = inv(W)
        for N in Ns
            # basis actually diagonalizes all Ns
            m = W * N * iW
            @test isdiag(m)
        end
    end

    tbl = Characters.CharacterTable(F, G, ccG)
    @test Characters.irreducible_characters(tbl) isa
          Vector{<:Characters.Character}

    chars = Characters.irreducible_characters(tbl)

    # checking the degrees
    @test sum(Int.(SymbolicWedderburn.degree.(chars)) .^ 2) == order(G)

    # orthogonality of characters over F
    @test [dot(χ, ψ) for χ in chars, ψ in chars] == I
end

function generictest_dixon_C(G, p = Characters.dixon_prime(G))
    F = Characters.FiniteFields.GF{p}
    ccG = Characters.conjugacy_classes(G)
    tbl = Characters.CharacterTable(F, G, ccG)
    chars_Fp = Characters.irreducible_characters(tbl)

    degrees = SymbolicWedderburn.degree.(chars_Fp)

    m = Characters._multiplicities(chars_Fp)
    for i in 1:size(m, 1)
        @test all(Int.(m[i, :, :]) .<= degrees[i])
    end

    @test Characters.complex_character_table(Rational{Int}, tbl) isa
          Characters.CharacterTable{<:Group,<:Cyclotomic{Rational{Int}}}

    tblCC = Characters.complex_character_table(Rational{Int}, tbl)
    chars_CC = Characters.irreducible_characters(tblCC)

    id = one(G)
    @test [ψ(id) for ψ in chars_CC] == degrees

    # orthogonality of characters over ℂ
    @test [dot(χ, ψ) for χ in chars_CC, ψ in chars_CC] == I

    @test sum(χ -> SymbolicWedderburn.degree(χ)^2, chars_CC) == order(G)
end

@testset "Dixon Algorithm" begin
    @testset "DixonPrimes" begin
        @test Characters.dixon_prime(20, 20) == 41

        @testset "random" begin
            for n in rand(2:1000, 10)
                F = Characters.Primes.factor(n)
                e = lcm(collect(keys(F)))
                p = Characters.dixon_prime(n, e)
                @test mod(p, e) == 1
                @test p > 2 * floor(sqrt(n))
            end
        end

        @testset "DixonPrimeGroups" begin
            G = PermGroup(perm"(1,2)", perm"(1,2,3,4)")
            ccG = Characters.conjugacy_classes(G)
            @test exponent(G) == 12
            @test exponent(ccG) == 12
            @test Characters.dixon_prime(G) == Characters.dixon_prime(ccG)
        end
    end

    @testset "Dixon over GF{p}" begin
        @testset "Example: Alt(4)" begin
            G = PermGroup([perm"(1,2,3)(4)", perm"(2,3,4)"])
            ccG = let S = gens(G)
                [
                    Orbit(one(G), S),
                    Orbit(G(perm"(2,3,4)"), S),
                    Orbit(G(perm"(2,4,3)"), S),
                    Orbit(G(perm"(1,2)(3,4)"), S),
                ] # the order is taken from GAP
            end

            let ccG = ccG, p = 29
                F = Characters.FiniteFields.GF{p}
                Ns = [Characters.CMMatrix(ccG, i) for i in 1:length(ccG)]
                @test_throws AssertionError Characters.common_esd(Ns, F)
            end

            let ccG = ccG, p = 31
                F = Characters.FiniteFields.GF{p}
                Ns = [Characters.CMMatrix(ccG, i) for i in 1:length(ccG)]
                esd = Characters.EigenSpaceDecomposition(F.(Ns[1]))
                esd = Characters.refine(esd, F.(Ns[2]))
                esd = Characters.refine(esd, F.(Ns[3]))
                @test isdiag(esd)

                @test isdiag(Characters.common_esd(Ns, F))
            end

            # generic tests
            generictest_dixon_Fp(G)

            # tests specific to G
            p = Characters.dixon_prime(ccG)
            F = Characters.FiniteFields.GF{p}

            tbl = Characters.CharacterTable(F, G)
            chars = Characters.irreducible_characters(tbl)

            @test sort(SymbolicWedderburn.degree.(chars)) == [1, 1, 1, 3]

            @test Set(Int.(values(χ)) for χ in chars) ==
                  Set([[3, 0, 0, 6], [1, 4, 2, 1], [1, 2, 4, 1], [1, 1, 1, 1]])
        end
    end

    @testset "Lifting characters to ℂ" begin
        @testset "Example: Alt(4)" begin
            G = PermGroup([perm"(1,2,3)(4)", perm"(2,3,4)"])
            ccG = let S = gens(G)
                [
                    Orbit(one(G), S),
                    Orbit(G(perm"(2,3,4)"), S),
                    Orbit(G(perm"(2,4,3)"), S),
                    Orbit(G(perm"(1,2)(3,4)"), S),
                ] # the order of cclasses is taken from GAP
            end

            @testset "PowerMap" begin
                pmG = Characters.PowerMap(ccG)
                @test size(pmG) == (4, 6)
                @test pmG[:, 0] == ones(Int, 4)
                @test pmG[:, 1] == 1:4
                @test pmG[:, 2] == [1, 3, 2, 1]
                @test pmG[:, 3] == [1, 1, 1, 4]
                @test pmG[:, 4] == [1, 2, 3, 1]
                @test pmG[:, 5] == [1, 3, 2, 4]
                @test all(isone.(first.(ccG) .^ 6))
            end

            # generic tests
            generictest_dixon_C(G)

            chars_C = Characters.irreducible_characters(G, ccG)
            E = Characters.Cyclotomics.E

            @test [collect(values(χ)) for χ in chars_C] == [
                E(3, 0) .* [3, 0, 0, -1],
                E(3, 0) .* [1, E(3, 2), E(3, 1), 1],
                E(3, 0) .* [1, E(3, 1), E(3, 2), 1],
                E(3, 0) .* [1, 1, 1, 1],
            ]
        end

        @testset "Example: Sym(4)" begin
            G = PermGroup([perm"(1,2)", perm"(1,2,3,4)"])
            ccG = let S = gens(G)
                ccG = [
                    Orbit(one(G), S),
                    Orbit(G(perm"(1,2)(4)"), S),
                    Orbit(G(perm"(1,2)(3,4)"), S),
                    Orbit(G(perm"(1,2,3)(4)"), S),
                    Orbit(G(perm"(1,2,3,4)"), S),
                ] # the order of cclasses is taken from GAP
            end
            generictest_dixon_Fp(G)
            generictest_dixon_C(G)

            chars = Characters.irreducible_characters(G, ccG)

            @test sort(SymbolicWedderburn.degree.(chars)) == [1, 1, 2, 3, 3]
            @test [collect(values(χ)) for χ in chars] == [
                [2, 0, 2, -1, 0],
                [3, 1, -1, 0, -1],
                [1, 1, 1, 1, 1],
                [1, -1, 1, 1, -1],
                [3, -1, -1, 0, 1],
            ]
        end

        @testset "example: C₅ ⋊ C₄" begin
            G = PermGroup(perm"(1,2,4,5,3)", perm"(2,5,3,4)")
            S = gens(G)

            ccG = [
                Orbit(one(G), S),
                Orbit(G(perm"(2,3)(4,5)"), S),
                Orbit(G(perm"(2,4,3,5)"), S),
                Orbit(G(perm"(2,5,3,4)"), S),
                Orbit(G(perm"(1,2,4,5,3)"), S),
            ]
            # the order of cclasses is taken from GAP
            generictest_dixon_Fp(G)
            generictest_dixon_C(G)

            chars = Characters.irreducible_characters(G, ccG)

            @test sort(SymbolicWedderburn.degree.(chars)) == [1, 1, 1, 1, 4]
            @test [collect(values(χ)) for χ in chars] == [
                E(4, 0) .* [4, 0, 0, 0, -1],
                E(4, 0) .* [1, 1, 1, 1, 1],
                E(4, 0) .* [1, 1, -1, -1, 1],
                E(4, 0) .* [1, -1, E(4), -E(4), 1],
                E(4, 0) .* [1, -1, -E(4), E(4), 1],
            ]
        end
    end

    @testset "Different base rings for characters_dixon" begin
        G = PermGroup(perm"(1,2,3,4)")
        @test eltype(Characters.irreducible_characters(Rational{Int}, G)) <:
              Characters.Character{<:Cyclotomic{Rational{Int}}}
        @test eltype(Characters.irreducible_characters(Rational{BigInt}, G)) <:
              Characters.Character{<:Cyclotomic{Rational{BigInt}}}
    end

    @time @testset "SmallPermGroups" begin
        for (ord, groups) in SmallPermGroups
            @testset "SmallGroup($ord, $n)" for (n, G) in enumerate(groups)
                @test Characters.irreducible_characters(G) isa
                      Vector{<:Characters.Character}
                generictest_dixon_Fp(G)
                generictest_dixon_C(G)
            end
        end
    end
end
