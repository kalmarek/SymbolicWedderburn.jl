function _test_gf(n::Int)
    
    F = collect(GF{n})
    
    @test characteristic(F[1]) == n
    @test int.(F) == collect(0:n-1)
    @test iszero(F[1])
    @test all(F .== F .+ F[1])
    @test all(F .* F[1] .== zeros(GF{n}, n))

    for f in F
        @test iszero(sum(f for _ in 1:characteristic(GF{n})))
        @test iszero(f*characteristic(GF{n}))
    end

    @test isone(F[2])
    @test all(F .== F .* F[2])
    @test iszero(F[2] + F[end])

    for f in F 
        if iszero(f)
            @test_throws DomainError inv(f)
        else
            @test isone(inv(f)*f)
        end
    end
end

@testset "PrimeFields" begin
    for p in [2, 3, 7]
        @testset "GF{$(p)}" begin
            _test_gf(p)
        end
    end
end

