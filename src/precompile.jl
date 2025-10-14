import PrecompileTools

struct _Word{T}
    alphabet::Vector{T}
    letters::Vector{Int}
end

function Base.:(==)(w::_Word, v::_Word)
    return w.alphabet == v.alphabet && w.letters == v.letters
end

function Base.hash(w::_Word, h::UInt = UInt(0))
    return hash(w.alphabet, hash(w.letters, hash(_Word, h)))
end

struct _OnLetters <: ByPermutations end

function action(::_OnLetters, p::AP.AbstractPermutation, w::_Word)
    return _Word(w.alphabet, [l^p for l in w.letters])
end

PrecompileTools.@setup_workload begin
    A = [:a, :b, :c]
    words = [_Word(A, [i]) for i in 1:3]
    for r in 2:4
        append!(
            words,
            [_Word(A, collect(w)) for w in Iterators.product(fill(1:3, r)...)],
        )
    end
    act = _OnLetters()
    PrecompileTools.@compile_workload begin
        G = PG.PermGroup(PG.perm"(1,2,3)", PG.perm"(1,2)")
        wd = WedderburnDecomposition(Rational{Int}, G, act, words, words)
        wdfl = WedderburnDecomposition(Float64, G, act, words, words)
    end
end
