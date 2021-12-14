struct CharacterTable{Gr,T,O} <: AbstractMatrix{T}
    group::Gr
    conjugacy_classes::Vector{O}
    inv_of::Vector{Int}
    pmap::PowerMap
    values::Matrix{T}
        # rows indexed by irreducible characters,
        # columns indexed by conjugacy classes
end

## Array Interface:
Base.size(chtbl::CharacterTable) = size(chtbl.values)
Base.getindex(chtbl::CharacterTable, i::Integer, j::Integer) = chtbl.values[i, j]

## Accessors
powermap(chtbl::CharacterTable) = chtbl.pmap
Base.parent(chtbl::CharacterTable) = chtbl.group

conjugacy_classes(chtbl::CharacterTable) = chtbl.conjugacy_classes
nconjugacy_classes(chtbl::CharacterTable) = size(chtbl.values, 2)
nirreps(chtbl::CharacterTable) = size(chtbl.values, 1)

## irreps
irreducible_characters(chtbl::CharacterTable) =
    [Character(chtbl, i) for i in 1:size(chtbl, 1)]
irreducible_characters(T::Type, chtbl::CharacterTable) =
    [Character{T}(chtbl, i) for i in 1:size(chtbl, 1)]
irreducible_characters(G::Group, cclasses = conjugacy_classes(G)) =
    irreducible_characters(Rational{Int}, G, cclasses)

function irreducible_characters(
    R::Type{<:Rational},
    G::Group,
    cclasses=conjugacy_classes(G)
)
    return irreducible_characters(CharacterTable(R, G, cclasses))
end

function trivial_character(chtbl::CharacterTable)
	# return Character(chtbl, findfirst(r->all(isone, r), eachrow(chtbl)))
	# can't use findfirst(f, eachrow(...)) on julia-1.6
	for i in 1:size(chtbl, 1)
		all(isone, @view(chtbl[i, :])) && return Character(chtbl, i)
	end
	# never hit, to keep compiler happy
	return Character(chtbl, 0)
end


## construcing tables

function CharacterTable(
    Fp::Type{<:FiniteFields.GF},
    G::Group,
    cclasses = conjugacy_classes(G),
)
    Ns = [CMMatrix(cclasses, i) for i in 1:length(cclasses)]
    esd = common_esd(Ns, Fp)
    @assert isdiag(esd)

    tbl = CharacterTable(G, cclasses, _inv_of(cclasses), PowerMap(cclasses), esd.basis)

    tbl = normalize!(tbl)
    return tbl
end

function CharacterTable(R::Type{<:Rational}, G::Group, cclasses=conjugacy_classes(G))
    Fp = FiniteFields.GF{dixon_prime(cclasses)}
    tblFp = CharacterTable(Fp, G, cclasses)
    return complex_character_table(R, tblFp)
end

function complex_character_table(
    R::Type{<:Rational},
    tblFp::CharacterTable{<:Group, <:FiniteFields.GF},
)
    charsFp = irreducible_characters(tblFp)
    mult_c = _multiplicities(charsFp)

    e = size(mult_c, 3) # the exponent

    C = Cyclotomics.Cyclotomic{R,Cyclotomics.SparseVector{R,Int}}
    values = Matrix{C}(undef, size(tblFp))

    Es = [E(e, k) for k in 0:e - 1]
    Threads.@threads for j in 1:size(tblFp, 2) # conjugacy_classes
        for i in 1:size(tblFp, 1) # characters
            # reduced_embedding may prevent overflow sometimes
            values[i, j] = Cyclotomics.reduced_embedding(
                sum(mult_c[i, j, k + 1] * Es[k + 1] for k in 0:e - 1)
            )
        end
    end

    return CharacterTable(
        parent(tblFp),
        conjugacy_classes(tblFp),
        tblFp.inv_of,
        powermap(tblFp),
        values,
    )
end

function _inv_of(cc::AbstractVector{<:AbstractOrbit})
    inv_of = zeros(Int, size(cc))
    for (i, c) in enumerate(cc)
        g = inv(first(c))
        inv_of[i] = something(findfirst(k -> g in k, cc), 0)
    end
    any(iszero, inv_of) && throw(
        ArgumentError(
            "Could not find the conjugacy class for inverse of $(first(cc[findfirst(iszero, inv_of)])).",
        ),
    )
    return inv_of
end

function normalize!(chtbl::CharacterTable{<:Group,<:FiniteFields.GF})
    id = one(parent(chtbl))
    for (i, χ) in enumerate(irreducible_characters(chtbl))
        k = χ(id)
        if !isone(k)
            chtbl.values[i, :] .*= inv(k)
        end
        # ⟨χ, χ⟩ = 1/d²

        deg = sqrt(inv(dot(χ, χ)))
        @debug "normalizing with" dot(χ, χ) χ(id) χ

        # normalizing χ
        chtbl.values[i, :] .*= deg
    end
    return chtbl
end

## "fancy" io

function Base.show(io::IO, ::MIME"text/plain", chtbl::CharacterTable)
    println(io, "Character table of ", parent(chtbl), " over $(eltype(chtbl))")

	if !(get(io, :limit, false)::Bool)
		screenheight = screenwidth = typemax(Int)
	else
		sz = displaysize(io)::Tuple{Int,Int}
		screenheight, screenwidth = sz[1] - 4, sz[2]
	end

    nirr = nirreps(chtbl)
    nccl = nconjugacy_classes(chtbl)
    pre_width = 1 + length(string(nirr))
    sep = " │ "

    if nirr > screenheight - 1
        rows_to_print = [1:screenheight - 2; nirr]
    else
        rows_to_print = [1:nirr;]
    end

    maxpossiblecols = div(screenwidth - pre_width - 3, 3)
    if nccl > maxpossiblecols
        cols_to_print = [1:maxpossiblecols - 1; nccl]
    else
	    cols_to_print = [1:nccl;]
    end

    A = Base.alignment(
        io,
        chtbl.values,
        rows_to_print,
        cols_to_print,
        screenwidth,
        screenwidth,
        2,
    )

	hellipsis = nccl > length(A) ? " …" : ""

    print(io, " "^(pre_width), "   ")
    Base.print_matrix_row(
        io,
        reshape(cols_to_print, (1, length(cols_to_print))),
        A,
        1,
        cols_to_print,
        " ",
    )

    println(io, hellipsis)
    println(io, "─"^pre_width, "─┬─", "─"^(sum(sum, A) + length(A) - 1), hellipsis)

    if nirr > screenheight - 1
        for i in 1:screenheight - 2
            print(io, rpad("χ$(FiniteFields.subscriptify(i))", pre_width), sep)
            Base.print_matrix_row(io, chtbl.values, A, i, cols_to_print, " ")
            println(io, hellipsis)
        end
        print(io, " ⋮", " "^(pre_width - 2), sep)
		Base.print_matrix_vdots(io, "⋮", A, "  ", 2, 1, false)
		println(io)

        print(io, rpad("χ$(FiniteFields.subscriptify(nirr))", pre_width), sep)
        Base.print_matrix_row(io, chtbl.values, A, nirr, cols_to_print, " ")
		print(io, hellipsis)
    else
        for i in 1:nirr
            print(io, rpad("χ$(FiniteFields.subscriptify(i))", pre_width), sep)
            Base.print_matrix_row(io, chtbl.values, A, i, cols_to_print, " ")
            i != nirr && println(io, hellipsis)
        end
    end
end
