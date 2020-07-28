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
        G = SymbolicWedderburn.AbstractAlgebra.SymmetricGroup(4)
        C = SymbolicWedderburn.PermutationGroups.conjugacy_classes(G)
        _ccmatrix(C)
        ccm = SymbolicWedderburn.CCMatrix(C, 2)
        

    end
    @testset "random" begin
        for i in 2:6
            G = SymbolicWedderburn.AbstractAlgebra.SymmetricGroup(i)
            _ccmatrix(SymbolicWedderburn.PermutationGroups.conjugacy_classes(G))
        end
    end

        @warn "Should actually be tested on subgroups of S_n"

end

