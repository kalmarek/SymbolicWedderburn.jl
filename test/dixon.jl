@testset "DixonPrime" begin
    @testset "DixonPrimeNumbers" begin
        @test SymbolicWedderburn.dixon_prime(20, 20) == 41
    end

    @testset "DixonPrimeGroups" begin
    end
end
@warn "This section needs more tests!"

