function _dixon_Fp(G, p = SymbolicWedderburn.dixon_prime(G))
    F = SymbolicWedderburn.FiniteFields.GF{p}
    ccG = conjugacy_classes(G)
    Ns = (SymbolicWedderburn.CCMatrix(ccG, i) for i = 1:length(ccG))
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
    @test sum(degree.(chars).^2) == order(G)

    # orthogonality of characters over F
    @test [dot(χ, ψ) for χ in chars, ψ in chars] == I
end

function _dixon_C(G, p=SymbolicWedderburn.dixon_prime(G))
    F = SymbolicWedderburn.FiniteFields.GF{p}
    ccG = conjugacy_classes(G)
    chars_Fp = SymbolicWedderburn.characters_dixon(ccG, F)

    degrees = degree.(chars_Fp)

    m = SymbolicWedderburn._multiplicities(chars_Fp)
    for i in 1:size(m,1)
        @test all(Int.(m[i, :, :]) .<= degrees[i])
    end

    @test SymbolicWedderburn.complex_characters(chars_Fp) isa Vector{<:SymbolicWedderburn.Character}
    @test length(SymbolicWedderburn.complex_characters(chars_Fp[1:1])) == 1

    chars_CC = SymbolicWedderburn.complex_characters(chars_Fp)
    id = one(first(first(ccG)))
    @test [ψ(id) for ψ in chars_CC] == degrees

    # orthogonality of characters over ℂ
    @test [dot(χ, ψ) for χ in chars_CC, ψ in chars_CC] == I
end

@testset "DixonPrime" begin
    @testset "DixonPrimes" begin
        @test SymbolicWedderburn.dixon_prime(20, 20) == 41

        @testset "random" for n in rand(2:1000, 10)
            F = SymbolicWedderburn.Primes.factor(n)
            e = lcm(collect(keys(F)))
            p = SymbolicWedderburn.dixon_prime(n, e)
            @test mod(p, e) == 1
            @test p > 2 * floor(sqrt(n))
        end

        @testset "DixonPrimeGroups" begin
            G = SymbolicWedderburn.AbstractAlgebra.SymmetricGroup(4)
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
                Ns = [SymbolicWedderburn.CCMatrix(ccG, i) for i = 1:length(ccG)]
                @test_throws AssertionError SymbolicWedderburn.common_esd(Ns, F)
            end

            let ccG = ccG, p = 31
                F = SymbolicWedderburn.FiniteFields.GF{p}
                Ns = [SymbolicWedderburn.CCMatrix(ccG, i) for i = 1:length(ccG)]
                esd = SymbolicWedderburn.EigenSpaceDecomposition(F.(Ns[1]))
                esd = SymbolicWedderburn.refine(esd, F.(Ns[2]))
                esd = SymbolicWedderburn.refine(esd, F.(Ns[3]))
                @test isdiag(esd)

                @test isdiag(SymbolicWedderburn.common_esd(Ns, F))
            end

            # generic tests
            _dixon_Fp(G)

            # tests specific to G
            p = SymbolicWedderburn.dixon_prime(ccG)
            F = SymbolicWedderburn.FiniteFields.GF{p}

            chars = SymbolicWedderburn.characters_dixon(ccG, F)

            @test sort(degree.(chars)) == [1, 1, 1, 3]

            @test Set(Int.(χ.vals) for χ in chars) == Set([
                [3, 0, 0, 6],
                [1, 4, 2, 1],
                [1, 2, 4, 1],
                [1, 1, 1, 1],
            ])

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
                ] # the order is taken from GAP
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
            _dixon_C(G)

            chars_C = SymbolicWedderburn.characters_dixon(G)
            E = SymbolicWedderburn.Cyclotomics.E

            @test Set([χ.vals for χ in chars_C]) == Set([
                E(3,0)*[3, 0, 0, -1],
                [E(3,0), E(3,1), E(3,2), E(3,0)],
                [E(3,0), E(3,2), E(3,1), E(3,0)],
                E(3,0)*[1, 1, 1,  1],
            ])



        end
    end
end
