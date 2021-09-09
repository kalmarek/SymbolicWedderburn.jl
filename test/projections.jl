function test_orthogonality(chars)
    projs = [SymbolicWedderburn.matrix_projection(χ) for χ in chars]

    res = Matrix{Bool}(undef, length(projs), length(projs))

    for j in axes(res, 2)
        for i in axes(res, 1)
            p, q = projs[i], projs[j]
            if p == q
                # p^2 == p
                res[i,j] = all(x -> isapprox(0.0, x; atol = 1e-12), p*p - p)
            else
                res[i,j] = all(x -> isapprox(0.0, x; atol = 1e-12), p*q)
            end
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
          PermutationGroups.degree(G)

    G = PermGroup([perm"(1,2,3,4)", perm"(1,2)"])
    chars = SymbolicWedderburn.irreducible_characters(G)
    @test test_orthogonality(chars)
    @test sum(first ∘ size, SymbolicWedderburn.image_basis.(chars)) ==
          PermutationGroups.degree(G)

    G = PermGroup([perm"(1,2,3)", perm"(2,3,4)"])
    chars = SymbolicWedderburn.irreducible_characters(G)
    @test test_orthogonality(chars)
    @test sum(first ∘ size, SymbolicWedderburn.image_basis.(chars)) ==
          PermutationGroups.degree(G)

    @time for ord = 2:16
        # @testset "SmallGroup($ord, $n)"
        for (n, G) in enumerate(SmallPermGroups[ord])
            chars = SymbolicWedderburn.irreducible_characters(Rational{Int}, G)

            @test test_orthogonality(chars)
            @test sum(
                first ∘ size,
                SymbolicWedderburn.image_basis.(chars),
            ) == PermutationGroups.degree(G)
        end
    end
end
