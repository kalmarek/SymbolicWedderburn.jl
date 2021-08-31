function test_orthogonality(chars)
    projs = [SymbolicWedderburn.matrix_projection(χ) for χ in chars]

    res = map(Iterators.product(projs, projs)) do (p, q)
        if p == q
            # p^2 == p
            all(x -> isapprox(0.0, x; atol = 1e-12), p*p - p)
        else
            res = p * q
            all(x -> isapprox(0.0, x; atol = 1e-12), p*q)
        end
    end
    return all(res)
end

@testset "Orthogonality of projections" begin
    G = PermGroup([perm"(1,2,3,4)"])
    chars = let irr = SymbolicWedderburn.irreducible_characters(G)
        if !all(isreal, irr)
            irr, _ = SymbolicWedderburn.affordable_real(irr)
        end
        SymbolicWedderburn.Character{Float64}.(irr)
    end
    @test test_orthogonality(chars)
    @test sum(first ∘ size, SymbolicWedderburn.image_basis.(chars)) ==
          degree(G)

    G = PermGroup([perm"(1,2,3,4)", perm"(1,2)"])
    chars = SymbolicWedderburn.irreducible_characters(G)
    @test test_orthogonality(chars)
    @test sum(first ∘ size, SymbolicWedderburn.image_basis.(chars)) ==
          degree(G)

    G = PermGroup([perm"(1,2,3)", perm"(2,3,4)"])
    chars = SymbolicWedderburn.irreducible_characters(G)
    @test test_orthogonality(chars)
    @test sum(first ∘ size, SymbolicWedderburn.image_basis.(chars)) ==
          degree(G)

    @time for ord = 2:16
        # @testset "SmallGroup($ord, $n)"
        for (n, G) in enumerate(SmallPermGroups[ord])
            chars = SymbolicWedderburn.irreducible_characters(Rational{Int}, G)

            @test test_orthogonality(chars)
            @test sum(
                first ∘ size,
                SymbolicWedderburn.image_basis.(chars),
            ) == degree(G)
        end
    end
end
