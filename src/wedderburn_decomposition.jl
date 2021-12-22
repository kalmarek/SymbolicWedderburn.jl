struct WedderburnDecomposition{B,iV,DS<:DirectSummand,Hom}
    basis::B
    invariants::iV
    Uπs::Vector{DS}
    hom::Hom
end

function WedderburnDecomposition(
    T::Type,
    G::Group,
    action::Action,
    basis_full,
    basis_half,
    S = Rational{Int};
    semisimple = false,
)
    tbl = CharacterTable(S, G)
    ehom = CachedExtensionHomomorphism(G, action, basis_half, precompute = true)

    Uπs = symmetry_adapted_basis(T, tbl, ehom; semisimple = semisimple)

    basis = StarAlgebras.Basis{UInt32}(basis_full)
    invariants = invariant_vectors(tbl, action, basis)

    return WedderburnDecomposition(basis, invariants, Uπs, ehom)
end

function Base.show(io::IO, wbdec::SymbolicWedderburn.WedderburnDecomposition)
    ds = direct_summands(wbdec)
    simple = all(issimple.(ds))
    dims = size.(ds, 1)
    norbs = length(invariant_vectors(wbdec))

    print(io, "Wedderburn Decomposition into $norbs orbits and $(length(ds))")
    all(simple) && print(io, " simple")
    println(io, " summands of sizes")
    return print(io, dims)
end

invariant_vectors(wbdec::WedderburnDecomposition) = wbdec.invariants
StarAlgebras.basis(wbdec::WedderburnDecomposition) = wbdec.basis
direct_summands(wbdec::WedderburnDecomposition) = wbdec.Uπs
Base.eltype(wbdec::WedderburnDecomposition) =
    eltype(eltype(direct_summands(wbdec)))

_tmps(wbdec::WedderburnDecomposition) =
    zeros.(eltype(wbdec), size.(direct_summands(wbdec)))

_fillfactor(M::AbstractMatrix) = count(!iszero, M) / length(M)
_fillfactor(M::AbstractSparseMatrix) = nnz(M) / length(M)

function diagonalize(
    M::AbstractMatrix,
    wbdec::WedderburnDecomposition,
    tmps = _tmps(wbdec),
)
    # return [degree(Uπ).*(Uπ*(M*Uπ')) for Uπ in summands(wbdec)]

    T = eltype(eltype(direct_summands(wbdec)))
    Mπs = [(d = size(Uπ, 1); zeros(T, d, d)) for Uπ in direct_summands(wbdec)]
    return diagonalize!(Mπs, M, direct_summands(wbdec), tmps)
end

function diagonalize!(
    Mπs,
    M::AbstractMatrix,
    wbdec::WedderburnDecomposition,
    tmps = _tmps(wbdec),
)
    return diagonalize!(Mπs, M, direct_summands(wbdec), tmps)
end

function diagonalize!(
    Mπs::AbstractVector{<:AbstractMatrix},
    M::AbstractMatrix,
    Uπs::AbstractVector{<:DirectSummand},
    tmps,
)
    @assert length(Mπs) == length(Uπs)

    for (π, Uπ) in enumerate(Uπs)
        imUπ = image_basis(Uπ) # Base.Matrix to allow BLAS paths below
        LinearAlgebra.mul!(tmps[π], imUπ, M)
        LinearAlgebra.mul!(Mπs[π], tmps[π], imUπ')
        zerotol!(Mπs[π], atol = eps(eltype(imUπ)) * max(size(imUπ)...))
        Mπs[π] .*= degree(Uπ)
    end
    return Mπs
end

function invariant_vectors(
    tbl::Characters.CharacterTable,
    act::Action,
    basis::StarAlgebras.Basis,
)
    trχ = Characters.Character{Rational{Int}}(Characters.trivial_character(tbl))

    ehom = ExtensionHomomorphism(act, basis)
    img = image_basis(ehom, trχ)

    # change the format of invariant_vectors to image_basis(ehom, trχ)
    return sparsevec.(eachrow(img))
end

function invariant_vectors(
    tbl::Characters.CharacterTable,
    act::ByPermutations,
    basis::StarAlgebras.Basis{T,I},
) where {T,I}
    G = parent(tbl)
    tovisit = trues(length(basis))
    invariant_vs = Vector{SparseVector{Rational{Int}}}()

    ordG = order(Int, G)
    elts = collect(G)
    orbit = zeros(I, ordG)

    for i in eachindex(basis)
        if tovisit[i]
            bi = basis[i]
            Threads.@threads for j in eachindex(elts)
                orbit[j] = basis[action(act, elts[j], bi)]
            end
            tovisit[orbit] .= false
            push!(
                invariant_vs,
                sparsevec(orbit, fill(1 // ordG, ordG), length(basis)),
            )
        end
    end
    return invariant_vs
end
