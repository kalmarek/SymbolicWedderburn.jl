import PrecompileTools

PrecompileTools.@setup_workload begin
    include(joinpath(@__DIR__, "..", "test", "free_words.jl"))

    struct OnLetters <: ByPermutations end
    function action(::OnLetters, p::AP.AbstractPermutation, w::Word)
        return Word(w.alphabet, [w.letters[i]^p for i in eachindex(w.letters)])
    end

    M = FreeWords([:a, :b, :c])
    words = collect(Iterators.take(M, nwords(M, 4)))
    act = OnLetters()

    PrecompileTools.@compile_workload begin
        G = PG.PermGroup(PG.perm"(1,2,3)", PG.perm"(1,2)")
        wd = WedderburnDecomposition(Rational{Int}, G, act, words, words)
        wdfl = WedderburnDecomposition(Float64, G, act, words, words)
    end
end
