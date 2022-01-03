struct WedderburnDecomposition{B,iV,DS<:DirectSummand,Hom}
    basis::B
    invariants::iV
    Uπs::Vector{DS}
    hom::Hom
end

"""
    WedderburnDecomposition([T::Type, ]G::Group, action::Action, basis_full, basis_half[, S; semisimple])
Compute `WedderburnDecomposition` related to `G` acting on `basis_half`, representing a form of Wedderburn-Artin decomposition.

This object is intended to be used to simplify problems of positive-semidefinite optimization.
# Arguments
 * `basis_full` corresponds to the basis (or index set) for the objective functional (the set of constraints);
 * `basis_half` corresponds to the basis (or index set) of the PSD constraint.

For description of the remaining arguments see [symmetry_adapted_basis](@ref).

Return a `wd::WedderburnDecomposition` object which defines:
 * `basis(wd)`: the original `basis_full`.
 * `invariant_vectors(wd)`: a basis for the subspace `⟨basis_full⟩^G` invariant under the action of `G`.
 * `direct_summands(wd)`: a vector of `DirectSummands` defining a map from `⟨basis_half⟩` to a direct (orthogonal) sum of subspaces
 * `diagonalize(M::AbstractMatrix, wd::WedderburnDecomposition)`: a map implementing the Wedderburn-Artin decomposition for matrices `M ∈ End(⟨basis_half⟩)` (i.e. with rows and columns indexed by the elements of `basis_half`).

See also: [symmetry_adapted_basis](@ref).
"""
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
        imUπ = convert(Matrix, image_basis(Uπ)) # Base.Matrix to allow BLAS paths below
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
    triv_χ = Characters.Character{Rational{Int}}(Characters.trivial_character(tbl))
    ehom =
        CachedExtensionHomomorphism(parent(tbl), act, basis, precompute = true)
    # ehom = ExtensionHomomorphism(act, basis)

    mpr = matrix_projection_irr(ehom, triv_χ)
    mpr, pivots = row_echelon_form!(mpr)
    img = mpr[1:length(pivots), :]

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
