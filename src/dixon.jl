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
    itr = iterate(Ns)
    @assert itr !== nothing
    esd = EigenSpaceDecomposition(F.(first(itr)))
    for N in Iterators.rest(Ns, last(itr))
        esd = refine(esd, F.(N))
        @debug N esd.eigspace_ptrs
        isdiag(esd) && return esd
    end
    return esd
end

characters_dixon(G::AbstractPermutationGroup) =
    characters_dixon(conjugacy_classes(G))

function characters_dixon(cclasses::AbstractVector{<:AbstractOrbit})
    p = dixon_prime(cclasses)
    chars_𝔽p = characters_dixon(cclasses, FiniteFields.GF{p})
    return complex_characters(chars_𝔽p)
end

function characters_dixon(
    cclasses::AbstractVector{<:AbstractOrbit},
    F::Type{<:FiniteFields.GF},
)
    Ns = [CMMatrix(cclasses, i) for i = 1:length(cclasses)]
    esd = common_esd(Ns, F)
    @assert isdiag(esd) "Class Matricies failed to diagonalize! $esd"
    inv_ccls = _inv_of(cclasses)
    return [
        normalize!(Character(vec(eigensubspace), inv_ccls, cclasses))
        for eigensubspace in esd
    ]
end

function _multiplicities(
    chars::AbstractVector{<:Character{F}},
    cclasses = conjugacy_classes(first(chars)),
) where {F<:FiniteFields.GF}

    e = Int(exponent(cclasses))
    ie = inv(F(e))
    ω = FiniteFields.rootofunity(F, e)


    multiplicities = zeros(Int, length(chars), length(cclasses), e)
    powermap = PowerMap(cclasses)
    for (i, χ) in enumerate(chars)
        for j = 1:length(cclasses), k = 0:e-1
            multiplicities[i, j, k+1] =
                Int(ie * sum(χ[powermap[j, l]] * ω^-(k * l) for l = 0:e-1))
        end
    end

    return multiplicities
end


function complex_characters(
    chars::AbstractVector{<:Character{F}},
) where {F<:FiniteFields.GF}

    cclasses = conjugacy_classes(first(chars))
    lccl = length(cclasses)
    mult_c = _multiplicities(chars, cclasses)
    e = size(mult_c, 3) # the exponent

    inv_of_cls = first(chars).inv_of

    C = Cyclotomics.Cyclotomic{Int,Cyclotomics.SparseVector{Int,Int}}
    # C = typeof(Cyclotomics.E(5))

    complex_chars = Vector{Character{C,eltype(cclasses)}}(undef, length(chars))


    for i = 1:length(complex_chars)
        complex_chars[i] = Character(
            [
                Cyclotomics.reduced_embedding(sum(
                    mult_c[i, j, k+1] * E(e, k) for k = 0:e-1
                )) for j = 1:lccl
            ],
            inv_of_cls,
            cclasses,
        )
    end

    return complex_chars
end
