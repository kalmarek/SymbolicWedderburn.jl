function _ccmatrix(C)
    r = length(C)
    for i = 1:r
        ccm = SymbolicWedderburn.CCMatrix(C, i)
        l = length(C[i])
        @test all([sum(ccm[:,j]) == l for j = 1:r])
    end
    ccm1 = SymbolicWedderburn.CCMatrix(C, 1)
    @test all([isone(ccm1[i,i]) for i = 1:r])
end

@testset "ConjClassMatrix" begin
    @testset "example" begin
        # S_4
        G = SymbolicWedderburn.AbstractAlgebra.SymmetricGroup(4)
        C = SymbolicWedderburn.PermutationGroups.conjugacy_classes(G)
        _ccmatrix(C)
        @test Matrix(SymbolicWedderburn.CCMatrix(C, 2))  == [0 1 0 0 0;
                                                             6 0 3 0 2;
                                                             0 4 0 4 0;
                                                             0 0 3 0 4;
                                                             0 1 0 2 0]
        @test Matrix(SymbolicWedderburn.CCMatrix(C, 3))  == [0 0 1 0 0;
                                                             0 4 0 4 0;
                                                             8 0 4 0 8;
                                                             0 4 0 4 0;
                                                             0 0 3 0 0]
        @test Matrix(SymbolicWedderburn.CCMatrix(C, 4))  == [0 0 0 1 0;
                                                             0 0 3 0 4;
                                                             0 4 0 4 0;
                                                             6 0 3 0 2;
                                                             0 2 0 1 0]
        @test Matrix(SymbolicWedderburn.CCMatrix(C, 5))  == [0 0 0 0 1;
                                                             0 1 0 2 0;
                                                             0 0 3 0 0;
                                                             0 2 0 1 0;
                                                             3 0 0 0 2]
        
        # Alt(4)
        a = SymbolicWedderburn.PermutationGroups.perm"(1,2,3)"
        b = SymbolicWedderburn.PermutationGroups.perm"(1,2,4)"
        G =  SymbolicWedderburn.PermutationGroups.PermGroup([a,b])
        C = SymbolicWedderburn.PermutationGroups.conjugacy_classes(G)
        _ccmatrix(C)
         @test_broken Matrix(SymbolicWedderburn.CCMatrix(C, 2)) == [0 1 0 0;
                                                             3 2 0 0;
                                                             0 0 3 0;
                                                             0 0 0 3]
         @test_broken Matrix(SymbolicWedderburn.CCMatrix(C, 3)) == [0 0 1 0;
                                                             0 0 3 0;
                                                             0 0 0 4;
                                                             4 4 0 0]
         @test_broken Matrix(SymbolicWedderburn.CCMatrix(C, 4)) == [0 0 0 1;
                                                             0 0 0 3;
                                                             4 4 0 0;
                                                             0 0 4 0]

        # I'd like to reproduce the example on p. 259.
            
    end
    @testset "random" begin
        for i in 2:6
            G = SymbolicWedderburn.AbstractAlgebra.SymmetricGroup(i)
            _ccmatrix(SymbolicWedderburn.PermutationGroups.conjugacy_classes(G))
        end
    end
end
@warn "Should actually be tested on subgroups of S_n"

