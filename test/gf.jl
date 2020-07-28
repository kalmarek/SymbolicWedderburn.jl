@testset "FiniteFields $p" for p in [2, 3, 7]

    @test GF{p}(0) isa GF{p}
    @test_throws AssertionError GF{-p}(0)
    @test_throws AssertionError GF{p^2}(0)

    @test int(zero(GF{p})) == 0
    @test int(one(GF{p})) == 1
    @test characteristic(GF{p}) == p
    @test characteristic(zero(GF{p})) == p

    @test GF{p}(0) == zero(GF{p})
    @test zero(GF{p}) != one(GF{p})
    @test iszero(GF{p}(0))
    @test !isone(GF{p}(0))
    @test isone(GF{p}(1))
    @test !isone(GF{p}(2))

    @test GF{p}(2p-1) == GF{p}(p-1) == GF{p}(-1)

    @test hash(GF{p}(1)) == hash(GF{p}(p+1)) == hash(one(GF{p}))
    @test hash(GF{p}(0)) != hash(GF{p}(1))

    p != 23 && @test_throws DomainError GF{p}(1) == GF{23}(1)

    @test -GF{p}(23) == GF{p}(-(23 % p))

    @test inv(GF{p}(23)) isa GF{p}
    @test isone(GF{p}(23)*inv(GF{p}(23)))
    @test_throws DomainError inv(GF{p}(p))

    @test all(GF{p}(i) + GF{p}(j) == GF{p}(i + j) for i in 0:p-1 for j in 0:p-1)
    @test all(GF{p}(i) - GF{p}(j) == GF{p}(i - j) for i in 0:p-1 for j in 0:p-1)
    @test all(GF{p}(i) * GF{p}(j) == GF{p}(i * j) for i in 0:p-1 for j in 0:p-1)

    @test all(iszero, (-GF{p}(i) + GF{p}(i) for i in 0:p-1))
    @test all(isone, (GF{p}(i)*inv(GF{p}(i)) for i in 1:p-1))
    @test all(GF{p}(i)^4 == GF{p}(i^4) for i in 0:p-1)
    @test all(GF{p}(23)/GF{p}(i) isa GF{p} for i in 1:p-1)

    @test all(GF{p}(i) + j == GF{p}(i + j) for i in 0:p-1 for j in 0:p-1)
    @test all(GF{p}(i) - j == GF{p}(i - j) for i in 0:p-1 for j in 0:p-1)
    @test all(GF{p}(i) * j == GF{p}(i * j) for i in 0:p-1 for j in 0:p-1)
    @test all(i/GF{p}(j) == GF{p}(i)/GF{p}(j) for i in 0:p-1 for j in 1:p-1)

    p!= 23 && @test_throws DomainError GF{p}(1) + GF{23}(2)
    p!= 23 && @test_throws DomainError GF{p}(1) - GF{23}(2)
    p!= 23 && @test_throws DomainError GF{p}(1) * GF{23}(2)
    p!= 23 && @test_throws DomainError GF{p}(1) / GF{23}(2)

    for i in 0:p
        x = GF{p}(i)
        try
            w = sqrt(x)
            @test w*w == x
            @test FiniteFields.issquare(x)
        catch er
            if er isa DomainError
                @test FiniteFields.legendresymbol(int(x), p) == -1
                @test !FiniteFields.issquare(x)
            else
                rethrow(er)
            end
        end
    end

    @test FiniteFields.generator(GF{p}) isa GF{p}
    g = FiniteFields.generator(GF{p})
    @test !iszero(g)
    @test sort(int.(g^i for i in 1:p-1)) == 1:p-1

    @test collect(GF{p}) isa Vector{GF{p}}
    F = collect(GF{p})
    @test int.(F) == collect(0:p-1)
    @test all(characteristic.(F) .== p)
    @test iszero(first(F))
    @test all(F .== F .+ F[1])
    @test all(F .* F[1] .== zeros(GF{p}, p))

    for f in F
        @test iszero(sum(f for _ in 1:characteristic(GF{p})))
        @test iszero(f*characteristic(GF{p}))
    end

    @test isone(F[2])
    @test all(F .== F .* F[2])

    @test string(one(GF{p})) == "1"*FiniteFields.subscriptify(p)
end
