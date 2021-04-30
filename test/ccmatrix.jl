function generic_tests_ccmatrix(C)
    r = length(C)
    for i = 1:r
        ccm = SymbolicWedderburn.CMMatrix(C, i)
        l = length(C[i])
        @test all([sum(ccm[:,j]) == l for j = 1:r])
    end
    ccm1 = SymbolicWedderburn.CMMatrix(C, 1)
    @test all([isone(ccm1[i,i]) for i = 1:r])
end

@testset "ConjClassMatrix" begin
    @testset "example: Symmetric group (4)" begin
        G = PermutationGroups.SymmetricGroup(4)

        @test conjugacy_classes(G) isa AbstractVector{<:AbstractOrbit}

        for cc in conjugacy_classes(G)
            @test all(permtype(g) == permtype(first(cc)) for g in cc)
        end

        C = conjugacy_classes(G)
        generic_tests_ccmatrix(C)
        @test SymbolicWedderburn.CMMatrix(C, 2)  == [0 1 0 0 0;
                                                     6 0 3 0 2;
                                                     0 4 0 4 0;
                                                     0 0 3 0 4;
                                                     0 1 0 2 0]
        @test SymbolicWedderburn.CMMatrix(C, 3)  == [0 0 1 0 0;
                                                     0 4 0 4 0;
                                                     8 0 4 0 8;
                                                     0 4 0 4 0;
                                                     0 0 3 0 0]
        @test SymbolicWedderburn.CMMatrix(C, 4)  == [0 0 0 1 0;
                                                     0 0 3 0 4;
                                                     0 4 0 4 0;
                                                     6 0 3 0 2;
                                                     0 2 0 1 0]
        @test SymbolicWedderburn.CMMatrix(C, 5)  == [0 0 0 0 1;
                                                     0 1 0 2 0;
                                                     0 0 3 0 0;
                                                     0 2 0 1 0;
                                                     3 0 0 0 2]
    end
    @testset "example: Alternating group (4)" begin
        a = perm"(1,2,3)(4)"
        b = perm"(1,2,4)"
        G = PermGroup([a,b])
        C = [
            Orbit(gens(G), one(G)),
            Orbit(gens(G), G(perm"(1,2)(3,4)")),
            Orbit(gens(G), G(perm"(1,2,3)(4)")),
            Orbit(gens(G), G(perm"(1,3,2)(4)")),
            ]

        generic_tests_ccmatrix(C)
        @test SymbolicWedderburn.CMMatrix(C, 1) == [1 0 0 0
                                                    0 1 0 0
                                                    0 0 1 0
                                                    0 0 0 1]
        @test SymbolicWedderburn.CMMatrix(C, 2) == [0 1 0 0
                                                    3 2 0 0
                                                    0 0 3 0
                                                    0 0 0 3]
        @test SymbolicWedderburn.CMMatrix(C, 3) == [0 0 1 0
                                                    0 0 3 0
                                                    0 0 0 4
                                                    4 4 0 0]
        @test SymbolicWedderburn.CMMatrix(C, 4) == [0 0 0 1
                                                    0 0 0 3
                                                    4 4 0 0
                                                    0 0 4 0]
    end

    @testset "example: C₅ ⋊ C₄" begin
        G = PermGroup([perm"(1,2,4,5,3)", perm"(2,5,3,4)"]);
        S = gens(G)

        ccG = [
            Orbit(S, one(G)),
            Orbit(S, G(perm"(2,3)(4,5)")),
            Orbit(S, G(perm"(2,4,3,5)")),
            Orbit(S, G(perm"(2,5,3,4)")),
            Orbit(S, G(perm"(1,2,4,5,3)")),
            ]
        # the order of cclasses is taken from GAP

        @assert sum(length, ccG) == order(G)

        @test SymbolicWedderburn.CMMatrix(ccG, 1) ==
        [ 1  0  0  0  0
          0  1  0  0  0
          0  0  1  0  0
          0  0  0  1  0
          0  0  0  0  1 ]
        @test SymbolicWedderburn.CMMatrix(ccG, 2) ==
        [ 0  1  0  0  0
          5  0  0  0  5
          0  0  0  5  0
          0  0  5  0  0
          0  4  0  0  0 ]
        @test SymbolicWedderburn.CMMatrix(ccG, 3) ==
        [ 0  0  1  0  0
          0  0  0  5  0
          0  5  0  0  0
          5  0  0  0  5
          0  0  4  0  0 ]
        @test SymbolicWedderburn.CMMatrix(ccG, 4) ==
        [ 0  0  0  1  0
          0  0  5  0  0
          5  0  0  0  5
          0  5  0  0  0
          0  0  0  4  0 ]
        @test SymbolicWedderburn.CMMatrix(ccG, 5) ==
        [ 0  0  0  0  1
          0  4  0  0  0
          0  0  4  0  0
          0  0  0  4  0
          4  0  0  0  3 ]

        generic_tests_ccmatrix(ccG)
        generic_tests_ccmatrix(conjugacy_classes(G)) # might be in different order
    end

    @testset "random subgroups of SymetricGroup(N)" begin
        for i in 2:6
            G = PermutationGroups.SymmetricGroup(i)
            for _ in 1:5
                PG = PermGroup(rand(G, 2))
                generic_tests_ccmatrix(SymbolicWedderburn.conjugacy_classes(PG))
            end
        end
    end
end
