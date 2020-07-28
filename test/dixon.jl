function _dixon_prime(n::Int)
    F = SymbolicWedderburn.Primes.factor(n)
    e = lcm(collect(keys(F)))
    p = SymbolicWedderburn.dixon_prime(n, e)
    @test mod(p, e) == 1
    @test p > 2*floor(sqrt(n))
end


@testset "DixonPrime" begin
    @testset "DixonPrimeNumbers" begin
        @testset "examples" begin
            @test SymbolicWedderburn.dixon_prime(20, 20) == 41
        end
        
        @testset "random" begin
            for i in abs.(rand(Int8, 10)).+2
                _dixon_prime(i)
            end
        end
    end

    @testset "DixonPrimeGroups" begin
        G = SymbolicWedderburn.AbstractAlgebra.SymmetricGroup(4)
        ccG = SymbolicWedderburn.PermutationGroups.conjugacy_classes(G)
        @test exponent(G) == 12
        @test exponent(ccG) == 12
        @test SymbolicWedderburn.dixon_prime(G) == SymbolicWedderburn.dixon_prime(ccG)
    end
end
@warn "This section needs more tests!"
