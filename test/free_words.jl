struct _Word{T}
    alphabet::Vector{T}
    letters::Vector{Int}

    function _Word(a::AbstractVector{T}, l::AbstractVector{<:Integer}) where {T}
        all(i -> 1 <= i <= length(a), l) ||
            throw(ArgumentError("Invalid word over alphabet $a: $w"))
        return new{T}(a, l)
    end
end
Base.length(w::_Word) = length(w.letters)
function Base.show(io::IO, w::_Word)
    if isone(w)
        print(io, "(id)")
    else
        join(io, w.alphabet[w.letters], "·")
    end
end

function Base.:(==)(w::_Word, v::_Word)
    return w.alphabet == v.alphabet && w.letters == v.letters
end
Base.hash(w::_Word, h::UInt) = hash(w.alphabet, hash(w.letters, hash(_Word, h)))

Base.one(w::_Word) = _Word(w.alphabet, Int[])
Base.isone(w::_Word) = length(w) == 0

function Base.:*(w::_Word, z::_Word)
    @assert w.alphabet == z.alphabet
    return _Word(w.alphabet, [w.letters; z.letters])
end

function SA.star(w::_Word)
    # star(:a) = :b
    # star(:b) = :a
    # star(:c) = :c

    star_d = Dict(1 => 2, 2 => 1)

    newletters = [get(star_d, l, l) for l in Iterators.reverse(w.letters)]
    return _Word(w.alphabet, newletters)
end

Base.isless(w::_Word, v::_Word) = w.letters < v.letters

struct _FreeWords{T}
    alphabet::Vector{T}
end

Base.one(fw::_FreeWords) = _Word(fw.alphabet, Int[])

Base.eltype(::Type{_FreeWords{T}}) where {T} = _Word{T}
Base.IteratorSize(::Type{<:_FreeWords}) = Base.IsInfinite()
function Base.iterate(aw::_FreeWords)
    w = _Word(aw.alphabet, Int[])
    stack = [w]
    return w, (stack, 1)
end

function Base.iterate(aw::_FreeWords, state)
    stack, l = state
    if l > length(aw.alphabet)
        popfirst!(stack)
        l = 1
    end
    w = first(stack)
    nw = _Word(aw.alphabet, [w.letters; l])
    push!(stack, nw)
    return nw, (stack, l + 1)
end

nwords(M::_FreeWords, maxl::Integer) = nwords(M, 0, maxl)
function nwords(M::_FreeWords, minl::Integer, maxl::Integer)
    maxl < minl && return zero(maxl)
    k = oftype(maxl, length(M.alphabet))
    return sum(k^i for i in minl:maxl)
end

allwords(M::_FreeWords, radius) = collect(Iterators.take(M, nwords(M, radius)))

struct _OnLetters <: ByPermutations end
function action(::_OnLetters, p::AP.AbstractPermutation, w::_Word)
    return _Word(w.alphabet, [w.letters[i]^p for i in eachindex(w.letters)])
end
