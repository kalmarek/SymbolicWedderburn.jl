Base.parent(A::StarAlgebra{<:Group}) = A.object

function StarAlgebras.AlgebraElement(
    χ::AbstractClassFunction,
    RG::StarAlgebra{<:Group},
)
    G = parent(RG)
    @assert G === parent(χ)
    b = basis(RG)

    dim = degree(χ)
    ord = order(Int, G)
    c = dim // ord

    T = Base._return_type(*, Tuple{typeof(c), eltype(χ)})

    x = AlgebraElement(zeros(T, length(b)), RG)

    for (v, cc) in zip(values(χ), conjugacy_classes(χ))
        for g in cc
            x[inv(g)] = c * v
        end
    end
    return x
end

function algebra_elt_from_support(support, RG::StarAlgebra)
    b = basis(RG)
    I = [b[s] for s in support]
    V = fill(1, length(support))
    return AlgebraElement(sparsevec(I, V, length(b)), RG)
end
struct CyclicSubgroups{Gr, GEl}
    group::Gr
    seen::Dict{Int,Set{GEl}}
    min_order::Int
    max_order::Int

    function CyclicSubgroups(G::Group; min_order=1, max_order=order(Int, G))
        seen = Dict{Int, Set{eltype(G)}}()
        return new{typeof(G), eltype(G)}(G, seen, min_order, max_order)
    end
end

Base.eltype(citr::CyclicSubgroups) = valtype(citr.seen)
Base.IteratorSize(::CyclicSubgroups) = Base.SizeUnknown()

function Base.iterate(citr::CyclicSubgroups)
    g, state = iterate(citr.group) # g is identity here
    @assert isone(g)
    if citr.min_order ≤ 1 ≤ citr.max_order
        citr.seen[ord] = Set([g])
        return Set([g]), state
    end
    return iterate(citr, state)
end

function Base.iterate(citr::CyclicSubgroups, state)
    k = iterate(citr.group, state)
    k === nothing && return nothing
    g, state = k
    ord = order(Int, g)
    if citr.min_order ≤ ord ≤ citr.max_order
        if !haskey(citr.seen, ord)
            citr.seen[ord] = Set([g])
            return Set(g^i for i in 1:ord), state
        else
            if any(g^i in citr.seen[ord] for i in 1:ord-1)
                return iterate(citr, state)
            else
                push!(citr.seen[ord], g)
                return Set(g^i for i in 1:ord), state
            end
        end
    end
    return iterate(citr, state)
end

function small_idempotents(
    RG::StarAlgebra{<:Group},
    subgroups = CyclicSubgroups(parent(RG), min_order = 2),
)
    return (algebra_elt_from_support(H, RG) // length(H) for H in subgroups)
end

@static if VERSION < v"1.3.0"
    function (χ::Character)(α::AlgebraElement{<:StarAlgebra{<:Group}})
        @assert parent(χ) === parent(parent(α))
        return sum(α(g) * χ(g) for g in supp(α))
    end
else
    function (χ::AbstractClassFunction)(α::AlgebraElement{<:StarAlgebra{<:Group}})
        @assert parent(χ) === parent(parent(α))
        return sum(α(g) * χ(g) for g in supp(α))
    end
end

function rank_one_projection(χ::Character, RG::StarAlgebra{<:Group};
    idems = small_idempotents(RG), # idems are AlgebraElements over Rational{Int}
    iters=3
)
    degree(χ) == 1 && return one(RG, Rational{Int}), true

    for µ in idems
        isone(χ(µ)) && µ^2 == µ && return µ, true
    end

    for n in 2:iters
        for elts in Iterators.product(ntuple(i -> idems, n)...)
            µ = *(elts...)
            isone(χ(µ)) && µ^2 == µ && return µ, true
        end
    end
    @debug "Could not find minimal projection for $χ"
    return one(RG, Rational{Int}), false
end

function minimal_projection_system(
    chars::AbstractVector{<:AbstractClassFunction},
    RG::StarAlgebra{<:Group}
)
    res = fetch.([@spawn_compat rank_one_projection(χ, RG) for χ in chars])

    r1p, simple = first.(res), last.(res) # rp1 are sparse storage

    mps = [isone(µ) ? AlgebraElement(χ, RG) : µ*AlgebraElement(χ, RG)
        for (µ,χ) in zip(r1p, chars)] # dense storage

    return mps, simple
end
