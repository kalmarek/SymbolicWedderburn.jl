struct Word{T}
    alphabet::Vector{T}
    letters::Vector{Int}

    function Word(a::AbstractVector{T}, l::AbstractVector{<:Integer}) where {T}
        all(i -> 1 <= i <= length(a), l) ||
            throw(ArgumentError("Invalid word over alphabet $a: $w"))
        return new{T}(a, l)
    end
end
Base.show(io::IO, w::Word) = join(io, w.alphabet[w.letters], "·")

function Base.:(==)(w::Word, v::Word)
    return w.alphabet == v.alphabet && w.letters == v.letters
end
function Base.hash(w::Word, h::UInt = UInt(0))
    return hash(w.alphabet, hash(w.letters, hash(Word, h)))
end

struct OnLetters <: SymbolicWedderburn.ByPermutations end
function SymbolicWedderburn.action(
    ::OnLetters,
    p::PermutationGroups.AbstractPermutation,
    w::Word,
)
    return Word(w.alphabet, [w.letters[i]^p for i in eachindex(w.letters)])
end

function allwords(A, radius)
    words = [Word(A, [i]) for i in 1:length(A)]
    for r in 2:radius
        append!(
            words,
            [Word(A, collect(w)) for w in Iterators.product(fill(1:3, r)...)],
        )
    end
    return words
end

@testset "Extending homomorphism" begin
    words = let A = [:a, :b, :c], radius = 4
        w = Word(A, [1, 2, 3, 2, 1])

        # (a·b·c·b·a)^(2,3) == a·c·b·c·a
        @test SymbolicWedderburn.action(OnLetters(), perm"(2,3)", w) ==
              Word(A, [1, 3, 2, 3, 1])

        words = allwords(A, radius)
    end

    G = PermGroup(perm"(1,2,3)", perm"(1,2)") # G acts on words permuting letters
    action = OnLetters()
    tbl = SymbolicWedderburn.CharacterTable(Rational{Int}, G)
    ehom = SymbolicWedderburn.CachedExtensionHomomorphism(
        Int32,
        G,
        action,
        words;
        precompute = true,
    )
    @test all(g ∈ keys(ehom.cache) for g in G) # we actually cached
    @test typeof(SymbolicWedderburn.induce(ehom, one(G))) == Perm{Int32}

    ehom = SymbolicWedderburn.CachedExtensionHomomorphism(
        G,
        action,
        words;
        precompute = true,
    )
    @test all(g ∈ keys(ehom.cache) for g in G) # we actually cached
    @test typeof(SymbolicWedderburn.induce(ehom, one(G))) == Perm{UInt16}

    ψ = SymbolicWedderburn.action_character(ehom, tbl)
    @test SymbolicWedderburn.constituents(ψ) == [40, 22, 18]
    irr = SymbolicWedderburn.irreducible_characters(tbl)
    multips = SymbolicWedderburn.constituents(ψ)
    @test dot(SymbolicWedderburn.degree.(irr), multips) == length(words)
    simple = isone.(SymbolicWedderburn.degree.(irr))
    @test simple == [false, true, true]

    inv_vec = SymbolicWedderburn.invariant_vectors(
        tbl,
        action,
        SymbolicWedderburn.basis(ehom),
    )
    @test length(inv_vec) == 22
    @test eltype(eltype(inv_vec)) == Rational{Int}

    @testset "semisimple decomposition" begin
        let i = 1
            χ, m, s = irr[i], multips[i], simple[i]
            b = SymbolicWedderburn.image_basis(ehom, χ)
            @test size(b, 1) == SymbolicWedderburn.degree(χ) * m == 80
        end

        let i = 2
            χ, m, s = irr[i], multips[i], simple[i]
            b = SymbolicWedderburn.image_basis(ehom, χ)
            @test size(b, 1) == SymbolicWedderburn.degree(χ) * m == 22
        end

        let i = 3
            χ, m, s = irr[i], multips[i], simple[i]
            b = SymbolicWedderburn.image_basis(ehom, χ)
            @test size(b, 1) == SymbolicWedderburn.degree(χ) * m == 18
        end

        @test symmetry_adapted_basis(G, action, words; semisimple = true) isa
              Vector{
            <:SymbolicWedderburn.DirectSummand{<:SymbolicWedderburn.Cyclotomic},
        }
        @test symmetry_adapted_basis(
            Rational{Int},
            G,
            action,
            words;
            semisimple = true,
        ) isa Vector{<:SymbolicWedderburn.DirectSummand{Rational{Int}}}
        @test symmetry_adapted_basis(
            Float64,
            G,
            action,
            words;
            semisimple = true,
        ) isa Vector{<:SymbolicWedderburn.DirectSummand{Float64}}

        sa_basis = symmetry_adapted_basis(G, action, words; semisimple = true)

        @test [convert(Matrix{Rational{Int}}, b) for b in sa_basis] isa Vector{Matrix{Rational{Int}}}
        @test [convert(Matrix{Float64}, b) for b in sa_basis] isa Vector{Matrix{Float64}}

        @test length(sa_basis) == 3
        @test multiplicity.(sa_basis) == [40, 22, 18]
        @test SymbolicWedderburn.degree.(sa_basis) == [2, 1, 1]
        @test size.(sa_basis, 1) ==
              multips .* SymbolicWedderburn.degree.(irr) ==
              [80, 22, 18]
        @test sum(first ∘ size, sa_basis) == length(words)
    end

    @testset "simple decomposition" begin
        RG = let G = G
            b = StarAlgebras.Basis{UInt16}(collect(G))
            StarAlgebra(G, b, (length(b), length(b)))
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

        @test symmetry_adapted_basis(G, action, words; semisimple = false) isa
              Vector{
            <:SymbolicWedderburn.DirectSummand{<:SymbolicWedderburn.Cyclotomic},
        }
        @test symmetry_adapted_basis(
            Rational{Int},
            G,
            action,
            words;
            semisimple = false,
        ) isa Vector{<:SymbolicWedderburn.DirectSummand{Rational{Int}}}
        @test symmetry_adapted_basis(
            Float64,
            G,
            action,
            words;
            semisimple = false,
        ) isa Vector{<:SymbolicWedderburn.DirectSummand{Float64}}

        sa_basis = symmetry_adapted_basis(G, action, words)

        @test length(sa_basis) == 3
        @test multiplicity.(sa_basis) == [40, 22, 18]
        @test SymbolicWedderburn.degree.(sa_basis) == [2, 1, 1]
        @test all(issimple, sa_basis)
        @test size.(sa_basis, 1) == multips == [40, 22, 18]
    end
end
