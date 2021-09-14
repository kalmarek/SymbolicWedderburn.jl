using StarAlgebras

@testset "affordable real degrees/dot" begin
    G = SmallPermGroups[10][2] # C₂⊕C₅
    tbl = SymbolicWedderburn.CharacterTable(Rational{Int}, G)
    chars = SymbolicWedderburn.irreducible_characters(tbl)

    @test all(isone ∘ SymbolicWedderburn.degree, chars)
    @test all(χ -> isone(dot(χ, χ)), chars)
    chars_fl = SymbolicWedderburn.Character{ComplexF64}.(chars)
    @test all(isone ∘ SymbolicWedderburn.degree, chars_fl)
    @test all(χ -> isapprox(dot(χ, χ), 1), chars_fl)

    charsR, _ = SymbolicWedderburn.affordable_real(chars)
    @test SymbolicWedderburn.degree.(charsR) ==
        [dot(χ, χ) for χ in charsR] ==
        [1,1,2,2,2,2]
    charsR_fl = SymbolicWedderburn.Character{Float64}.(charsR)
    @test all([1,1,2,2,2,2] .==
        SymbolicWedderburn.degree.(charsR_fl) .≈
        [dot(χ, χ) for χ in charsR_fl]
        )
end

@testset "Symmetry adapted basis" begin
    @testset "step by step" begin
        G = PermGroup(perm"(1,2)", perm"(1,2,3)")
        irr = SymbolicWedderburn.irreducible_characters(G)
        @test irr isa AbstractVector{<:SymbolicWedderburn.Character{<:Cyclotomic}}
        @test SymbolicWedderburn.degree.(irr) == [2,1,1]

        RG = let G = G
            b = StarAlgebras.Basis{UInt16}(collect(G))
            StarAlgebra(G, b, (length(b), length(b)))
        end

        µ, s = SymbolicWedderburn.minimal_rank_projection(irr[1], RG)
        @test µ isa AlgebraElement
        @test s

        @test SymbolicWedderburn.Character{Rational{Int}}(irr[1]) isa
        SymbolicWedderburn.Character{Rational{Int}}
        @test collect(values(SymbolicWedderburn.Character{Float64}(irr[1]))) isa Vector{Float64}

        mps, simple = SymbolicWedderburn.minimal_projection_system(irr, RG)
        @test all(simple)

        @test rank(float.(SymbolicWedderburn.matrix_projection(irr[1]))) == 2
        @test rank(float.(SymbolicWedderburn.matrix_projection(mps[1]))) == 1

        sa_basis_ssimple = symmetry_adapted_basis(Rational{Int}, G, Rational{Int}, semisimple=true)

        @test issimple.(sa_basis_ssimple) == [false, true]
        @test rank.(convert.(Matrix, sa_basis_ssimple)) == [2, 1]
        @test dot(
            multiplicity.(sa_basis_ssimple),
            SymbolicWedderburn.degree.(sa_basis_ssimple)
        ) == PermutationGroups.degree(G)

        sa_basis = symmetry_adapted_basis(Rational{Int}, G, Rational{Int}, semisimple=false)
        @test all(issimple.(sa_basis))
        @test rank.(convert.(Matrix{Float64}, sa_basis)) == [1,1]
    end

    C₄ = PermGroup([perm"(1,2,3,4)"])
    A₄ = PermGroup([perm"(1,2,3)", perm"(2,3,4)"])
    S₄ = PermGroup([perm"(1,2,3,4)", perm"(1,2)"])


    sa_basis = symmetry_adapted_basis(C₄, semisimple=true)
    @test sum(first ∘ size, sa_basis) == PermutationGroups.degree(C₄)
    sa_basis = symmetry_adapted_basis(A₄, semisimple=true)
    @test sum(first ∘ size, sa_basis) == PermutationGroups.degree(A₄)
    sa_basis = symmetry_adapted_basis(S₄, semisimple=true)
    @test sum(first ∘ size, sa_basis) == PermutationGroups.degree(S₄)

    sa_basis = symmetry_adapted_basis(ComplexF64, C₄, semisimple=true)
    @test sum(first ∘ size, sa_basis) == PermutationGroups.degree(C₄)
    sa_basis = symmetry_adapted_basis(ComplexF64, A₄, semisimple=true)
    @test sum(first ∘ size, sa_basis) == PermutationGroups.degree(A₄)
    sa_basis = symmetry_adapted_basis(ComplexF64, S₄, semisimple=true)
    @test sum(first ∘ size, sa_basis) == PermutationGroups.degree(S₄)

    sa_basis = symmetry_adapted_basis(Float64, C₄, semisimple=true)
    @test sum(first ∘ size, sa_basis) == PermutationGroups.degree(C₄)
    sa_basis = symmetry_adapted_basis(Float64, A₄, semisimple=true)
    @test sum(first ∘ size, sa_basis) == PermutationGroups.degree(A₄)
    sa_basis = symmetry_adapted_basis(Float64, S₄, semisimple=true)
    @test sum(first ∘ size, sa_basis) == PermutationGroups.degree(S₄)
end

@testset "Symmetry adapted basis: small groups" begin
    @time @testset "Rational" begin
        for ord = 2:30
            @testset "SmallGroup($ord, $n)" for (n, G) in enumerate(SmallPermGroups[ord])
                @testset "complex" begin
                    S = Rational{Int}

                    sa_basis = symmetry_adapted_basis(G, S, semisimple=true)

                    @test dot(SymbolicWedderburn.degree.(sa_basis), multiplicity.(sa_basis)) ==
                        PermutationGroups.degree(G)

                    @test sum(first ∘ size, sa_basis) == PermutationGroups.degree(G)

                    S = if (ord, n) in ((21,1),)
                        Rational{BigInt}
                    else
                        Rational{Int}
                    end
                    sa_basis = symmetry_adapted_basis(G, S, semisimple=false)
                    for b in sa_basis
                        if issimple(b)
                            @test multiplicity(b) == size(b, 1)
                        else
                            @test multiplicity(b)*SymbolicWedderburn.degree(b) == size(b, 1)
                        end
                    end
                end

                @testset "real" begin
                    S = Rational{Int}

                    sa_basisR = symmetry_adapted_basis(Float64, G, S, semisimple=true)
                    @test dot(SymbolicWedderburn.degree.(sa_basisR), multiplicity.(sa_basisR)) ==
                        PermutationGroups.degree(G)
                    @test sum(first ∘ size, sa_basisR) == PermutationGroups.degree(G)

                    sa_basisR = symmetry_adapted_basis(Float64, G, S, semisimple=false)
                    for b in sa_basisR
                        if issimple(b)
                            @test multiplicity(b) == size(b, 1) || 2*multiplicity(b) == size(b, 1)
                            # the first condiditon doesn't hold for realified characters;
                        else
                            @test multiplicity(b)*SymbolicWedderburn.degree(b) == size(b, 1)
                        end
                    end
                end
            end
        end
    end

    @time @testset "Float64/Complex64" begin
        for ord = 2:30
            @testset "SmallGroup($ord, $n)" for (n, G) in enumerate(SmallPermGroups[ord])

                T = Float64
                sa_basis = symmetry_adapted_basis(T, G, semisimple=true)
                @test sa_basis isa Vector{<:AbstractMatrix{T}}

                @test all(sa_basis) do b
                    rank(b) == size(b, 1)
                end

                @test sum(first ∘ size, sa_basis) == PermutationGroups.degree(G)

                basisC = symmetry_adapted_basis(Complex{T}, G, semisimple=true)
                @test sum(first ∘ size, basisC) == PermutationGroups.degree(G)
            end
        end
    end
end
