@testset "Symmetry adapted basis" begin
    G = PermGroup([perm"(1,2,3,4)"])
    basis = symmetry_adapted_basis(G)
    @test sum(first ∘ size ∘ SymbolicWedderburn.basis, basis) == degree(G)

    G = PermGroup([perm"(1,2,3,4)", perm"(1,2)"])
    basis = symmetry_adapted_basis(ComplexF64, G)
    @test sum(first ∘ size ∘ SymbolicWedderburn.basis, basis) == degree(G)

    G = PermGroup([perm"(1,2,3)", perm"(2,3,4)"])
    basis = symmetry_adapted_basis(ComplexF64, G)
    @test sum(first ∘ size ∘ SymbolicWedderburn.basis, basis) == degree(G)

    G = PermGroup([perm"(1,2,3,4)"])
    basis = symmetry_adapted_basis(Float64, G)
    @test sum(first ∘ size ∘ SymbolicWedderburn.basis, basis) == degree(G)

    G = PermGroup([perm"(1,2,3,4)", perm"(1,2)"])
    basis = symmetry_adapted_basis(Float64, G)
    @test sum(first ∘ size ∘ SymbolicWedderburn.basis, basis) == degree(G)

    G = PermGroup([perm"(1,2,3)", perm"(2,3,4)"])
    basis = symmetry_adapted_basis(G)
    @test sum(first ∘ size ∘ SymbolicWedderburn.basis, basis) == degree(G)
end


@testset "Symmetry adapted basis: small groups" begin
    @time @testset "Rational" begin
        for ord = 2:30
            @testset "SmallGroup($ord, $n)" for (n, G) in enumerate(SmallPermGroups[ord])

                @testset "complex" begin
                    S = if (ord, n) in ((21,1),)
                        Rational{BigInt}
                    else
                        Rational{Int}
                    end

                    basis = symmetry_adapted_basis(G, S)

                    @test dot(SymbolicWedderburn.degree.(basis), SymbolicWedderburn.multiplicity.(basis)) == degree(G)

                    @test sum(first ∘ size ∘ SymbolicWedderburn.basis, basis) == degree(G)
                end

                @testset "real" begin
                    S = if (ord, n) in ((17,1), (19,1), (21,1), (23,1), (24,3), (25,1), (29,1))
                        Rational{BigInt}
                    else
                        Rational{Int}
                    end

                    chars = SymbolicWedderburn.affordable_real!(
                        SymbolicWedderburn.characters_dixon(S, G)
                    )

                    basisR = symmetry_adapted_basis(chars)
                    @test dot(SymbolicWedderburn.degree.(basisR), SymbolicWedderburn.multiplicity.(basisR)) == degree(G)

                    @test sum(first ∘ size ∘ SymbolicWedderburn.basis, basisR) == degree(G)
                end
            end
        end
    end

    @time @testset "Float64/Complex64" begin
        for ord = 2:30
            @testset "SmallGroup($ord, $n)" for (n, G) in enumerate(SmallPermGroups[ord])

                T = Float64
                basis = symmetry_adapted_basis(T, G)

                @test all(basis) do b
                    B = SymbolicWedderburn.basis(b)
                    rank(B) == size(B, 1)
                end

                @test sum(first ∘ size ∘ SymbolicWedderburn.basis, basis) == degree(G)

                basisC = symmetry_adapted_basis(Complex{T}, G)
                @test sum(first ∘ size ∘ SymbolicWedderburn.basis, basisC) == degree(G)
            end
        end
    end
end
