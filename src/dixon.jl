Base.exponent(G::AbstractPermutationGroup) = exponent(conjugacy_classes(G))
Base.exponent(cclasses::AbstractVector) = lcm(order.(first.(cclasses)))
dixon_prime(G::AbstractPermutationGroup) = dixon_prime(order(G), exponent(G))

function dixon_prime(cclasses::AbstractVector)
    ordG = sum(length, cclasses)
    m = exponent(cclasses)
    return dixon_prime(ordG, m)
end

function dixon_prime(ordG::Integer, exponent::Integer)
    p = 2 * floor(Int, sqrt(ordG))
    while true
        p = nextprime(p + 1)
        isone(p % exponent) && break # we need -1 to be in the field
    end
    return p
end

function common_esd(Ns, F::Type{<:FiniteFields.GF})
    @assert !isempty(Ns)
    esd = EigenSpaceDecomposition(F.(first(Ns)))
    for N in Iterators.rest(Ns, 2)
        esd = refine(esd, F.(N))
        @debug N esd.eigspace_ptrs
        isdiag(esd) && return esd
    end
    return esd
end

function characters_dixon(
    cclasses::AbstractVector{<:AbstractOrbit},
    F::Type{<:FiniteFields.GF},
)
    Ns = (CCMatrix(cclasses, i) for i = 1:length(cclasses))
    esd = common_esd(Ns, F)
    inv_ccls = _inv_of(cclasses)
    return [
        normalize!(Character(vec(eigensubspace), inv_ccls, cclasses))
        for eigensubspace in esd
    ]
end
