include("free_words.jl")

struct OnLetters <: SymbolicWedderburn.ByPermutations end
function SymbolicWedderburn.action(
    ::OnLetters,
    p::AP.AbstractPermutation,
    w::Word,
)
    return Word(w.alphabet, [w.letters[i]^p for i in eachindex(w.letters)])
end

struct OnLettersSigned <: SymbolicWedderburn.BySignedPermutations end
function SymbolicWedderburn.action(
    ::OnLettersSigned,
    p::AP.AbstractPermutation,
    w::Word,
)
    return (Word(w.alphabet, [w.letters[i]^p for i in eachindex(w.letters)]), 1)
end

@testset "Extending homomorphism" begin
    fb_words = let A = [:a, :b, :c], radius = 4
        w = Word(A, [1, 2, 3, 2, 1])

        # (a·b·c·b·a)^(2,3) == a·c·b·c·a
        @test SymbolicWedderburn.action(OnLetters(), PG.perm"(2,3)", w) ==
              Word(A, [1, 3, 2, 3, 1])

        StarAlgebras.FixedBasis(
            allwords(FreeWords(A), radius),
            StarAlgebras.DiracMStructure(*),
        )
    end

    G = PG.PermGroup(PG.perm"(1,2,3)", PG.perm"(1,2)") # G acts on words permuting letters

    @test SymbolicWedderburn.check_group_action(G, OnLetters(), fb_words)
    @test SymbolicWedderburn.check_group_action(G, OnLettersSigned(), fb_words)

    action = OnLetters()
    ehom = SymbolicWedderburn.ExtensionHomomorphism(action, fb_words)
    @test typeof(SymbolicWedderburn.induce(ehom, one(G))) == PG.Perm{UInt32}

    let T = UInt16, fb_words = fb_words
        l = length(fb_words)
        fb_words = StarAlgebras.FixedBasis(
            collect(fb_words),
            StarAlgebras.DiracMStructure(*),
            T.((l, l)),
        )
        @test SymbolicWedderburn.check_group_action(G, OnLetters(), fb_words)
        @test SymbolicWedderburn.check_group_action(
            G,
            OnLettersSigned(),
            fb_words,
        )

        action = OnLetters()
        ehom = SymbolicWedderburn.ExtensionHomomorphism(action, fb_words)
        @test typeof(SymbolicWedderburn.induce(ehom, one(G))) == PG.Perm{T}
    end

    schrhom = SymbolicWedderburn.SchreierExtensionHomomorphism(
        G,
        action,
        fb_words;
        memoize = false,
    )

    @test all(
        SymbolicWedderburn.induce(ehom, g) ==
        SymbolicWedderburn.induce(schrhom, g) for g in G
    )
    @test length(schrhom.cache) == ngens(G) + 1

    schrhom = SymbolicWedderburn.SchreierExtensionHomomorphism(
        G,
        action,
        fb_words;
        memoize = true,
    )

    @test all(
        SymbolicWedderburn.induce(ehom, g) ==
        SymbolicWedderburn.induce(schrhom, g) for g in G
    )
    @test length(schrhom.cache) == order(Int, G)

    tbl = SymbolicWedderburn.CharacterTable(Rational{Int}, G)
    ψ = SymbolicWedderburn.action_character(ehom, tbl)
    @test SymbolicWedderburn.multiplicities(ψ) == [23, 18, 40]
    ψ = SymbolicWedderburn.action_character(schrhom, tbl)
    @test SymbolicWedderburn.multiplicities(ψ) == [23, 18, 40]

    irr = SymbolicWedderburn.irreducible_characters(tbl)
    multips = SymbolicWedderburn.multiplicities(ψ)
    @test dot(SymbolicWedderburn.degree.(irr), multips) == length(basis(ehom))
    simple = isone.(SymbolicWedderburn.degree.(irr))
    @test simple == [true, true, false]

    inv_vec = SymbolicWedderburn.invariant_vectors(
        tbl,
        action,
        SymbolicWedderburn.basis(ehom),
    )
    @test length(inv_vec) == 23
    @test eltype(eltype(inv_vec)) == Rational{Int}

    @testset "semisimple decomposition" begin
        let i = 1
            χ, m, s = irr[i], multips[i], simple[i]
            b = SymbolicWedderburn.image_basis(schrhom, χ)
            @test size(b, 1) == SymbolicWedderburn.degree(χ) * m == 23
        end

        let i = 2
            χ, m, s = irr[i], multips[i], simple[i]
            b = SymbolicWedderburn.image_basis(schrhom, χ)
            @test size(b, 1) == SymbolicWedderburn.degree(χ) * m == 18
        end

        let i = 3
            χ, m, s = irr[i], multips[i], simple[i]
            b = SymbolicWedderburn.image_basis(schrhom, χ)
            @test size(b, 1) == SymbolicWedderburn.degree(χ) * m == 80
        end

        @test symmetry_adapted_basis(G, action, fb_words; semisimple = true) isa
              Vector{
            <:SymbolicWedderburn.DirectSummand{<:SymbolicWedderburn.Cyclotomic},
        }
        @test symmetry_adapted_basis(
            Rational{Int},
            G,
            action,
            fb_words;
            semisimple = true,
        ) isa Vector{<:SymbolicWedderburn.DirectSummand{Rational{Int}}}
        @test symmetry_adapted_basis(
            Float64,
            G,
            action,
            fb_words;
            semisimple = true,
        ) isa Vector{<:SymbolicWedderburn.DirectSummand{Float64}}

        sa_basis =
            symmetry_adapted_basis(G, action, fb_words; semisimple = true)

        @test [convert(Matrix{Rational{Int}}, b) for b in sa_basis] isa Vector{Matrix{Rational{Int}}}
        @test [convert(Matrix{Float64}, b) for b in sa_basis] isa Vector{Matrix{Float64}}

        @test length(sa_basis) == 3
        @test multiplicity.(sa_basis) == [23, 18, 40]
        @test SymbolicWedderburn.degree.(sa_basis) == [1, 1, 2]
        @test size.(sa_basis, 1) ==
              multips .* SymbolicWedderburn.degree.(irr) ==
              [23, 18, 80]
        @test sum(first ∘ size, sa_basis) == length(fb_words)
    end

    @testset "simple decomposition" begin
        RG = let G = G
            v = collect(G)
            l = convert(UInt16, length(v))
            b = StarAlgebras.FixedBasis(
                v,
                StarAlgebras.DiracMStructure(*),
                (l, l),
            )
            StarAlgebra(G, b)
        end

        let χ = irr[1], m = multips[1]
            (a, fl) = SymbolicWedderburn.minimal_rank_projection(χ, RG)
            @test fl == isone(χ(a))

            µ = AlgebraElement(χ, RG) * a
            mpr = SymbolicWedderburn.image_basis(ehom, µ)
            @test mpr isa AbstractMatrix{eltype(µ)}
            @test size(mpr, 1) == m

            µR = AlgebraElement{Rational{Int}}(µ)
            mpr = SymbolicWedderburn.image_basis(ehom, µR)
            @test mpr isa AbstractMatrix{eltype(µR)}
            @test size(mpr, 1) == m

            µFl = AlgebraElement{Float64}(µ)
            mpr = SymbolicWedderburn.image_basis(ehom, µFl)
            @test mpr isa AbstractMatrix{eltype(µFl)}
            @test size(mpr, 1) == m
        end

        let χ = irr[2], m = multips[2]
            (a, fl) = SymbolicWedderburn.minimal_rank_projection(χ, RG)
            @test fl == isone(χ(a))

            µ = AlgebraElement(χ, RG) * a
            mpr = SymbolicWedderburn.image_basis(ehom, µ)
            @test mpr isa AbstractMatrix{eltype(µ)}
            @test size(mpr, 1) == m
        end

        let χ = irr[3], m = multips[3]
            (a, fl) = SymbolicWedderburn.minimal_rank_projection(χ, RG)
            @test fl == isone(χ(a))

            µ = AlgebraElement(χ, RG) * a
            mpr = SymbolicWedderburn.image_basis(ehom, µ)
            @test mpr isa AbstractMatrix{eltype(µ)}
            @test size(mpr, 1) == m
        end

        @test symmetry_adapted_basis(
            G,
            action,
            fb_words;
            semisimple = false,
        ) isa
              Vector{
            <:SymbolicWedderburn.DirectSummand{<:SymbolicWedderburn.Cyclotomic},
        }
        @test symmetry_adapted_basis(
            Rational{Int},
            G,
            action,
            fb_words;
            semisimple = false,
        ) isa Vector{<:SymbolicWedderburn.DirectSummand{Rational{Int}}}
        @test symmetry_adapted_basis(
            Float64,
            G,
            action,
            fb_words;
            semisimple = false,
        ) isa Vector{<:SymbolicWedderburn.DirectSummand{Float64}}

        sa_basis = symmetry_adapted_basis(G, action, fb_words)

        @test length(sa_basis) == 3
        @test multiplicity.(sa_basis) == [23, 18, 40]
        @test SymbolicWedderburn.degree.(sa_basis) == [1, 1, 2]
        @test all(issimple, sa_basis)
        @test size.(sa_basis, 1) == multips == [23, 18, 40]
    end
end
