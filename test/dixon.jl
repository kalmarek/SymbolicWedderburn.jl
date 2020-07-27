@testset "DixonPrime" begin
    @testset "DixonPrimeNumbers" begin
        G = SymbolicWedderburn.AbstractAlgebra.SymmetricGroup(4)
        ccG = SymbolicWedderburn.PermutationGroups.conjugacy_classes(G)
        @test exponent(G) == 12
        @test exponent(ccG) == 12

        @test SymbolicWedderburn.dixon_prime(G) == SymbolicWedderburn.dixon_prime(ccG)
        @test SymbolicWedderburn.dixon_prime(20, 20) == 41
    end

    @testset "DixonPrimeGroups" begin
    end
end
@warn "This section needs more tests!"
