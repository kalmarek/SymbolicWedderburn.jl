struct Word{T}
    alphabet::Vector{T}
    letters::Vector{Int}

    function Word(a::AbstractVector{T}, l::AbstractVector{<:Integer}) where T
        all(i->1<=i<=length(a), l) || throw(ArgumentError("Invalid word over alphabet $a: $w"))
        return new{T}(a,l)
    end
end
Base.show(io::IO, w::Word) = join(io, w.alphabet[w.letters], "·")

Base.:(==)(w::Word, v::Word) = w.alphabet == v.alphabet && w.letters == v.letters
Base.hash(w::Word, h::UInt=UInt(0)) = hash(w.alphabet, hash(w.letters, hash(Word, h)))

struct OnLetters <: SymbolicWedderburn.ByPermutations end
function SymbolicWedderburn.action(a::OnLetters, p::PermutationGroups.AbstractPerm, w::Word)
    return Word(w.alphabet, [l^p for l in w.letters])
end

my_action(w::Word, g) = SymbolicWedderburn.action(OnLetters(), g, w)

@testset "Extending homomorphism" begin
    words = let A = [:a, :b, :c]
        w = Word(A, [1,2,3,2,1])

        # (a·b·c·b·a)^(2,3) == a·c·b·c·a
        @test my_action(w, perm"(2,3)") == Word(A, [1,3,2,3,1])

        words = [Word(A, [1]), Word(A, [2]), Word(A, [3])]
        for r in 2:4
            append!(
                words,
                [Word(A, collect(w)) for w in Iterators.product(fill(1:3, r)...)]
            )
        end
        words
    end

    G = PermGroup(perm"(1,2,3)",perm"(1,2)") # G acts on words permuting letters
    sa_basis = symmetry_adapted_basis(G, words, OnLetters())

    @test sa_basis isa Vector{<:SymbolicWedderburn.SemisimpleSummand{<:Matrix{<:SymbolicWedderburn.Cyclotomic}}}
    @test [convert(Matrix{Int}, b) for b in sa_basis] isa Vector{Matrix{Int}}
    @test [convert(Matrix{Float64}, b) for b in sa_basis] isa Vector{Matrix{Float64}}
    @test sum(first ∘ size ∘ SymbolicWedderburn.basis, sa_basis) == length(words)
end
