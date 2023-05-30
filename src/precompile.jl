import PrecompileTools

PrecompileTools.@setup_workload begin
    struct Word{T}
        alphabet::Vector{T}
        letters::Vector{Int}

        function Word(
            a::AbstractVector{T},
            l::AbstractVector{<:Integer},
        ) where {T}
            all(i -> 1 <= i <= length(a), l) ||
                throw(ArgumentError("Invalid word over alphabet $a: $l"))
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

    struct OnLetters <: ByPermutations end
    function action(::OnLetters, p::PermutationGroups.AbstractPerm, w::Word)
        return Word(w.alphabet, [w.letters[i]^p for i in eachindex(w.letters)])
    end

    function allwords(A, radius)
        words = [Word(A, [i]) for i in 1:length(A)]
        for r in 2:radius
            append!(
                words,
                [
                    Word(A, collect(w)) for
                    w in Iterators.product(fill(1:3, r)...)
                ],
            )
        end
        return words
    end

    words = allwords([:a, :b, :c], 4)
    act = OnLetters()

    PrecompileTools.@compile_workload begin
        G = PermGroup(perm"(1,2,3)", perm"(1,2)")
        wd = WedderburnDecomposition(Rational{Int}, G, act, words, words)
        wdfl = WedderburnDecomposition(Float64, G, act, words, words)
    end
end
