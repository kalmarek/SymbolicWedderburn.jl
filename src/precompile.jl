import PrecompileTools

include(joinpath(@__DIR__, "..", "test", "free_words.jl"))

PrecompileTools.@setup_workload begin
    M = _FreeWords([:a, :b, :c])
    words = collect(Iterators.take(M, nwords(M, 4)))
    act = _OnLetters()

    PrecompileTools.@compile_workload begin
        G = PG.PermGroup(PG.perm"(1,2,3)", PG.perm"(1,2)")
        wd = WedderburnDecomposition(Rational{Int}, G, act, words, words)
        wdfl = WedderburnDecomposition(Float64, G, act, words, words)
    end
end
