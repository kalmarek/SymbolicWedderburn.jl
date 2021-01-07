function generictest_dixon_Fp(G, p = SymbolicWedderburn.dixon_prime(G))
    F = SymbolicWedderburn.FiniteFields.GF{p}
    ccG = conjugacy_classes(G)
    Ns = [SymbolicWedderburn.CMMatrix(ccG, i) for i = 1:length(ccG)]
    @test isdiag(SymbolicWedderburn.common_esd(Ns, F))
    esd = SymbolicWedderburn.common_esd(Ns, F)

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

    @test SymbolicWedderburn.characters_dixon(ccG, F) isa
          Vector{<:SymbolicWedderburn.Character}

    chars = SymbolicWedderburn.characters_dixon(ccG, F)

    # checking the degrees
    @test sum(degree.(chars) .^ 2) == order(G)

    # orthogonality of characters over F
    @test [dot(χ, ψ) for χ in chars, ψ in chars] == I
end

function generictest_dixon_C(G, p = SymbolicWedderburn.dixon_prime(G))
    F = SymbolicWedderburn.FiniteFields.GF{p}
    ccG = conjugacy_classes(G)
    chars_Fp = SymbolicWedderburn.characters_dixon(ccG, F)

    degrees = degree.(chars_Fp)

    m = SymbolicWedderburn._multiplicities(chars_Fp)
    for i = 1:size(m, 1)
        @test all(Int.(m[i, :, :]) .<= degrees[i])
    end

    @test SymbolicWedderburn.complex_characters(chars_Fp) isa
          Vector{<:SymbolicWedderburn.Character}
    @test length(SymbolicWedderburn.complex_characters(chars_Fp[1:1])) == 1

    chars_CC = SymbolicWedderburn.complex_characters(chars_Fp)
    id = one(first(first(ccG)))
    @test [ψ(id) for ψ in chars_CC] == degrees

    # orthogonality of characters over ℂ
    @test [dot(χ, ψ) for χ in chars_CC, ψ in chars_CC] == I
end

@testset "Dixon Algorithm" begin
    @testset "DixonPrimes" begin
        @test SymbolicWedderburn.dixon_prime(20, 20) == 41

        @testset "random" begin
            for n in rand(2:1000, 10)
                F = SymbolicWedderburn.Primes.factor(n)
                e = lcm(collect(keys(F)))
                p = SymbolicWedderburn.dixon_prime(n, e)
                @test mod(p, e) == 1
                @test p > 2 * floor(sqrt(n))
            end
        end

        @testset "DixonPrimeGroups" begin
            G = PermutationGroups.SymmetricGroup(4)
            ccG = conjugacy_classes(G)
            @test exponent(G) == 12
            @test exponent(ccG) == 12
            @test SymbolicWedderburn.dixon_prime(G) ==
                  SymbolicWedderburn.dixon_prime(ccG)
        end
    end
    @testset "Dixon over GF{p}" begin
        @testset "Example: Alt(4)" begin
            G = PermGroup([perm"(1,2,3)(4)", perm"(2,3,4)"])
            ccG = let S = gens(G)
                [
                    Orbit(S, one(G)),
                    Orbit(S, perm"(2,3,4)"),
                    Orbit(S, perm"(2,4,3)"),
                    Orbit(S, perm"(1,2)(3,4)"),
                ] # the order is taken from GAP
            end

            let ccG = ccG, p = 29
                F = SymbolicWedderburn.FiniteFields.GF{p}
                Ns = [SymbolicWedderburn.CMMatrix(ccG, i) for i = 1:length(ccG)]
                @test_throws AssertionError SymbolicWedderburn.common_esd(Ns, F)
            end

            let ccG = ccG, p = 31
                F = SymbolicWedderburn.FiniteFields.GF{p}
                Ns = [SymbolicWedderburn.CMMatrix(ccG, i) for i = 1:length(ccG)]
                esd = SymbolicWedderburn.EigenSpaceDecomposition(F.(Ns[1]))
                esd = SymbolicWedderburn.refine(esd, F.(Ns[2]))
                esd = SymbolicWedderburn.refine(esd, F.(Ns[3]))
                @test isdiag(esd)

                @test isdiag(SymbolicWedderburn.common_esd(Ns, F))
            end

            # generic tests
            generictest_dixon_Fp(G)

            # tests specific to G
            p = SymbolicWedderburn.dixon_prime(ccG)
            F = SymbolicWedderburn.FiniteFields.GF{p}

            chars = SymbolicWedderburn.characters_dixon(ccG, F)

            @test sort(degree.(chars)) == [1, 1, 1, 3]

            @test Set(Int.(χ.vals) for χ in chars) ==
                  Set([[3, 0, 0, 6], [1, 4, 2, 1], [1, 2, 4, 1], [1, 1, 1, 1]])

        end
    end

    @testset "Lifting characters to ℂ" begin
        @testset "Example: Alt(4)" begin
            G = PermGroup([perm"(1,2,3)(4)", perm"(2,3,4)"])
            ccG = let S = gens(G)
                [
                    Orbit(S, one(G)),
                    Orbit(S, perm"(2,3,4)"),
                    Orbit(S, perm"(2,4,3)"),
                    Orbit(S, perm"(1,2)(3,4)"),
                ] # the order of cclasses is taken from GAP
            end

            @testset "PowerMap" begin
                pmG = SymbolicWedderburn.PowerMap(ccG)
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

            chars_C = SymbolicWedderburn.characters_dixon(G)
            E = SymbolicWedderburn.Cyclotomics.E

            @test [χ.vals for χ in chars_C] == [
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
                    Orbit(S, one(G)),
                    Orbit(S, perm"(1,2)(4)"),
                    Orbit(S, perm"(1,2)(3,4)"),
                    Orbit(S, perm"(1,2,3)(4)"),
                    Orbit(S, perm"(1,2,3,4)"),
                ] # the order of cclasses is taken from GAP
            end
            generictest_dixon_Fp(G)
            generictest_dixon_C(G)

            chars = SymbolicWedderburn.characters_dixon(ccG)

            @test sort(degree.(chars)) == [1, 1, 2, 3, 3]
            @test [χ.vals for χ in chars] == [
                [2, 0, 2, -1, 0],
                [3, 1, -1, 0, -1],
                [1, 1, 1, 1, 1],
                [1, -1, 1, 1, -1],
                [3, -1, -1, 0, 1],
            ]
        end

        @testset "example: C₅ ⋊ C₄" begin
            S = [perm"(1,2,4,5,3)", perm"(2,5,3,4)"]
            G = PermGroup(S)

            ccG = [
                Orbit(S, one(G)),
                Orbit(S, perm"(2,3)(4,5)"),
                Orbit(S, perm"(2,4,3,5)"),
                Orbit(S, perm"(2,5,3,4)"),
                Orbit(S, perm"(1,2,4,5,3)"),
            ]
            # the order of cclasses is taken from GAP
            generictest_dixon_Fp(G)
            generictest_dixon_C(G)

            chars = SymbolicWedderburn.characters_dixon(ccG)

            @test sort(degree.(chars)) == [1, 1, 1, 1, 4]
            @test [χ.vals for χ in chars] == [
                E(4, 0) .* [4, 0, 0, 0, -1],
                E(4, 0) .* [1, 1, 1, 1, 1],
                E(4, 0) .* [1, 1, -1, -1, 1],
                E(4, 0) .* [1, -1, E(4), -E(4), 1],
                E(4, 0) .* [1, -1, -E(4), E(4), 1],
            ]
        end
    end

    @testset "SmallPermGroups" begin
        include("smallgroups.jl")
        for (ord, groups) in SmallPermGroups
            @testset "SmallGroup($ord, $n)" for (n, G) in enumerate(groups)
                @test SymbolicWedderburn.characters_dixon(G) isa
                      Vector{<:SymbolicWedderburn.Character}
                generictest_dixon_Fp(G)
                generictest_dixon_C(G)
            end
        end
    end
end
#=
@testset "Characters io" begin
    G = PermGroup([perm"(2,3)(4,5)"])
    chars = SymbolicWedderburn.characters_dixon(G)

    @test sprint(show, chars[1]) == "SymbolicWedderburn.Character: [1, 1]"
    @test sprint(show, chars[2]) == "SymbolicWedderburn.Character: [1, -1]"

    @test sprint(show, MIME"text/plain"(), chars[1]) ==
        "SymbolicWedderburn.Character over Cyclotomic{Int64,SparseArrays.SparseVector{Int64,Int64}}\n()^G         → \t 1*E(1)^0\n(2,3)(4,5)^G → \t 1*E(1)^0"
    @test sprint(show, MIME"text/plain"(), chars[2]) ==
        "SymbolicWedderburn.Character over Cyclotomic{Int64,SparseArrays.SparseVector{Int64,Int64}}\n()^G         → \t 1*E(1)^0\n(2,3)(4,5)^G → \t-1*E(1)^0"
end
=#
