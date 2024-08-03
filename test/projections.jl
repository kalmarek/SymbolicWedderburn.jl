function test_orthogonality(chars)
    projs = [SymbolicWedderburn.matrix_projection(χ) for χ in chars]

    res = Matrix{Bool}(undef, length(projs), length(projs))

    for j in axes(res, 2)
        for i in axes(res, 1)
            p, q = projs[i], projs[j]
            if p == q
                # p^2 == p
                res[i, j] = all(x -> isapprox(0.0, x; atol = 1e-12), p * p - p)
            else
                res[i, j] = all(x -> isapprox(0.0, x; atol = 1e-12), p * q)
            end
        end
    end
    return all(res)
end

@testset "Orthogonality of projections" begin
    G = PermGroup([perm"(1,2,3,4)"])
    chars = let irr = SymbolicWedderburn.irreducible_characters(G)
        if !all(isreal, irr)
            irr, _ = SymbolicWedderburn.affordable_real(irr,fill(1, length(irr)))
        end
        SymbolicWedderburn.Character{Float64}.(irr)
    end
    @test test_orthogonality(chars)
    @test sum(first ∘ size, SymbolicWedderburn.image_basis.(chars)) ==
          AP.degree(G)

    G = PermGroup([perm"(1,2,3,4)", perm"(1,2)"])
    chars = SymbolicWedderburn.irreducible_characters(G)
    @test test_orthogonality(chars)
    @test sum(first ∘ size, SymbolicWedderburn.image_basis.(chars)) ==
          AP.degree(G)

    G = PermGroup([perm"(1,2,3)", perm"(2,3,4)"])
    chars = SymbolicWedderburn.irreducible_characters(G)
    @test test_orthogonality(chars)
    @test sum(first ∘ size, SymbolicWedderburn.image_basis.(chars)) ==
          AP.degree(G)

    for ord in 2:16
        for (n, G) in enumerate(SmallPermGroups[ord])
            @testset "SmallGroup($ord, $n)" begin
                chars =
                    SymbolicWedderburn.irreducible_characters(Rational{Int}, G)

                @test test_orthogonality(chars)

                @test sum(
                    first ∘ size ∘ SymbolicWedderburn.image_basis,
                    chars,
                ) == AP.degree(G)

                mps, _ = SymbolicWedderburn.minimal_projection_system(
                    chars,
                    SymbolicWedderburn._group_algebra(G),
                )
                @test sum(first ∘ size ∘ SymbolicWedderburn.image_basis, mps) ≤
                      AP.degree(G)
            end
        end
    end
end

@testset "affordable real all cases" begin
    G = SmallPermGroups[24][3] # SL(2,3)
    tbl = SymbolicWedderburn.CharacterTable(Rational{Int}, G)
    chars = SymbolicWedderburn.irreducible_characters(tbl)

    @test SymbolicWedderburn.degree.(chars)==[1, 1, 1, 2, 2, 2, 3]
    @test all(χ -> isone(dot(χ, χ)), chars)
    chars_fl = SymbolicWedderburn.Character{ComplexF64}.(chars)
    @test SymbolicWedderburn.degree.(chars_fl)==[1, 1, 1, 2, 2, 2, 3]
    @test all(χ -> isapprox(dot(χ, χ), 1), chars_fl)

  

    charsR, _ = SymbolicWedderburn.affordable_real(chars,fill(2,length(chars)))
    @test SymbolicWedderburn.degree.(charsR) == [1, 2, 4, 4, 3]
    @test [dot(χ, χ) for χ in charsR] == [1, 2, 2, 4, 1]
          
    charsR_fl = SymbolicWedderburn.Character{Float64}.(charsR)
    @test all( [1, 2, 4, 4, 3] .== SymbolicWedderburn.degree.(charsR_fl),)
    @test all( [1, 2, 2, 4, 1].≈[dot(χ, χ) for χ in charsR_fl],)
end
