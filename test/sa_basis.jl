@testset "affordable real degrees/dot" begin
    G = SmallPermGroups[10][2] # C₂⊕C₅
    tbl = SymbolicWedderburn.CharacterTable(Rational{Int}, G)
    chars = SymbolicWedderburn.irreducible_characters(tbl)

    @test all(isone ∘ SymbolicWedderburn.degree, chars)
    @test all(χ -> isone(dot(χ, χ)), chars)
    chars_fl = SymbolicWedderburn.Character{ComplexF64}.(chars)
    @test all(isone ∘ SymbolicWedderburn.degree, chars_fl)
    @test all(χ -> isapprox(dot(χ, χ), 1), chars_fl)

    charsR, _ = SymbolicWedderburn.affordable_real(chars,fill(1, length(chars)))
    @test SymbolicWedderburn.degree.(charsR) ==
          [dot(χ, χ) for χ in charsR] ==
          [1, 1, 2, 2, 2, 2]
    charsR_fl = SymbolicWedderburn.Character{Float64}.(charsR)
    @test all(
        [1, 1, 2, 2, 2, 2] .==
        SymbolicWedderburn.degree.(charsR_fl) .≈
        [dot(χ, χ) for χ in charsR_fl],
    )
end

@testset "Symmetry adapted basis" begin
    @testset "step by step" begin
        G = PermGroup(perm"(1,2)", perm"(1,2,3)")
        irr = SymbolicWedderburn.irreducible_characters(G)
        @test irr isa
              AbstractVector{<:SymbolicWedderburn.Character{<:Cyclotomic}}
        @test SymbolicWedderburn.degree.(irr) == [1, 1, 2]

        RG = let G = G
            l = order(UInt16, G)
            b = SA.FixedBasis(collect(G), SA.DiracMStructure(*), (l, l))
            StarAlgebra(G, b)
        end

        µ, s = SymbolicWedderburn.minimal_rank_projection(irr[1], RG)
        @test µ isa AlgebraElement
        @test isone(s)

        @test SymbolicWedderburn.Character{Rational{Int}}(irr[1]) isa
              SymbolicWedderburn.Character{Rational{Int}}
        @test collect(values(SymbolicWedderburn.Character{Float64}(irr[1]))) isa
              Vector{Float64}

        mps, ranks = SymbolicWedderburn.minimal_projection_system(irr, RG)
        @test all(isone, ranks)

        @test rank(float.(SymbolicWedderburn.matrix_projection(irr[3]))) == 2
        @test rank(float.(SymbolicWedderburn.matrix_representation(mps[3]))) ==
              1

        sa_basis_ssimple = symmetry_adapted_basis(
            Rational{Int},
            G,
            Rational{Int};
            semisimple = true,
        )

        @test issimple.(sa_basis_ssimple) == [true, false]
        @test rank.(convert.(Matrix, sa_basis_ssimple)) == [1, 2]
        @test dot(
            multiplicity.(sa_basis_ssimple),
            SymbolicWedderburn.degree.(sa_basis_ssimple),
        ) == AP.degree(G)

        sa_basis = symmetry_adapted_basis(
            Rational{Int},
            G,
            Rational{Int};
            semisimple = false,
        )
        @test all(issimple.(sa_basis))
        @test rank.(convert.(Matrix{Float64}, sa_basis)) == [1, 1]
    end

    C₄ = PermGroup([perm"(1,2,3,4)"])
    A₄ = PermGroup([perm"(1,2,3)", perm"(2,3,4)"])
    S₄ = PermGroup([perm"(1,2,3,4)", perm"(1,2)"])

    sa_basis = symmetry_adapted_basis(C₄; semisimple = true)
    @test sum(first ∘ size, sa_basis) == AP.degree(C₄)
    sa_basis = symmetry_adapted_basis(A₄; semisimple = true)
    @test sum(first ∘ size, sa_basis) == AP.degree(A₄)
    sa_basis = symmetry_adapted_basis(S₄; semisimple = true)
    @test sum(first ∘ size, sa_basis) == AP.degree(S₄)

    sa_basis = symmetry_adapted_basis(ComplexF64, C₄; semisimple = true)
    @test sum(first ∘ size, sa_basis) == AP.degree(C₄)
    sa_basis = symmetry_adapted_basis(ComplexF64, A₄; semisimple = true)
    @test sum(first ∘ size, sa_basis) == AP.degree(A₄)
    sa_basis = symmetry_adapted_basis(ComplexF64, S₄; semisimple = true)
    @test sum(first ∘ size, sa_basis) == AP.degree(S₄)

    sa_basis = symmetry_adapted_basis(Float64, C₄; semisimple = true)
    @test sum(first ∘ size, sa_basis) == AP.degree(C₄)
    sa_basis = symmetry_adapted_basis(Float64, A₄; semisimple = true)
    @test sum(first ∘ size, sa_basis) == AP.degree(A₄)
    sa_basis = symmetry_adapted_basis(Float64, S₄; semisimple = true)
    @test sum(first ∘ size, sa_basis) == AP.degree(S₄)
end

@testset "Symmetry adapted basis: small groups" begin
    @time @testset "Rational" begin
        for ord in 2:30
            @testset "SmallGroup($ord, $n)" for (n, G) in
                                                enumerate(SmallPermGroups[ord])
                @testset "complex" begin
                    S = Rational{Int}

                    sa_basis = symmetry_adapted_basis(G, S; semisimple = true)

                    @test dot(
                        SymbolicWedderburn.degree.(sa_basis),
                        multiplicity.(sa_basis),
                    ) == AP.degree(G)

                    @test sum(first ∘ size, sa_basis) ==
                          AP.degree(G)

                    S = if (ord, n) in ((26, 1),)
                        Rational{BigInt}
                    else
                        Rational{Int}
                    end
                    sa_basis = symmetry_adapted_basis(G, S; semisimple = false)
                    for b in sa_basis
                        if issimple(b)
                            @test multiplicity(b) == size(b, 1) ||
                                  2 * multiplicity(b) == size(b, 1)
                        else
                            @test multiplicity(b) *
                                  SymbolicWedderburn.degree(b) == size(b, 1)
                        end
                    end
                end

                @testset "real" begin
                    S = Rational{Int}

                    sa_basisR =
                        symmetry_adapted_basis(Float64, G, S; semisimple = true)
                    @test dot(
                        SymbolicWedderburn.degree.(sa_basisR),
                        multiplicity.(sa_basisR),
                    ) == AP.degree(G)
                    @test sum(first ∘ size, sa_basisR) ==
                          AP.degree(G)

                    sa_basisR = symmetry_adapted_basis(
                        Float64,
                        G,
                        S;
                        semisimple = false,
                    )
                    for b in sa_basisR
                        if issimple(b)
                            @test multiplicity(b) == size(b, 1) ||
                                  2 * multiplicity(b) == size(b, 1)
                            # the first condiditon doesn't hold for realified characters;
                        else
                            rk, res = divrem(size(b, 1), multiplicity(b))
                            @test res == 0
                            @test 1 ≤ rk ≤ SymbolicWedderburn.degree(b)
                        end
                    end
                end
            end
        end
    end

    @time @testset "Float64/Complex64" begin
        for ord in 2:30
            @testset "SmallGroup($ord, $n)" for (n, G) in
                                                enumerate(SmallPermGroups[ord])
                T = Float64
                sa_basis = symmetry_adapted_basis(T, G; semisimple = true)
                @test sa_basis isa Vector{<:AbstractMatrix{T}}

                @test all(sa_basis) do b
                    return rank(b) == size(b, 1)
                end

                @test sum(first ∘ size, sa_basis) == AP.degree(G)

                basisC =
                    symmetry_adapted_basis(Complex{T}, G; semisimple = true)
                @test sum(first ∘ size, basisC) == AP.degree(G)
            end
        end
    end
end
