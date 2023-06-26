struct WedderburnDecomposition{B,iV,DS<:DirectSummand,Hom}
    basis::B
    invariants::iV
    Uπs::Vector{DS}
    hom::Hom
end

"""
    WedderburnDecomposition([T::Type, ]G::Group, action::Action, basis_full, basis_half[, S; semisimple])
Compute `WedderburnDecomposition` related to `G` acting on `basis_half`, representing a form of Wedderburn-Artin decomposition.

This object is intended to be used to simplify problems of positive-semidefinite
optimization.
# Arguments
 * `basis_full` corresponds to the basis (or the index set) for the objective
   functional (the set of constraints);
 * `basis_half` corresponds to the basis (or the index set) for the PSD constraint.

For description of the remaining arguments see [symmetry_adapted_basis](@ref).

Return a `wd::WedderburnDecomposition` object which defines:
 * `basis(wd)`: the original `basis_full`.
 * `invariant_vectors(wd)`: a basis for the subspace `⟨basis_full⟩^G` invariant
    under the action of `G`.
 * `direct_summands(wd)`: a vector of `DirectSummands` defining a map from `⟨basis_half⟩`
    to a direct (orthogonal) sum of subspaces (the diagonalization map).
 * `diagonalize(M::AbstractMatrix, wd::WedderburnDecomposition)`: a map
    implementing the Wedderburn-Artin decomposition for matrices `M ∈ End(⟨basis_half⟩)`
    (i.e. with rows and columns indexed by the elements of `basis_half`).

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
    try_harder = false,
)
    tbl = CharacterTable(S, G)
    ehom = CachedExtensionHomomorphism(G, action, basis_half; precompute = true)
    check_group_action(G, ehom; full_check = false)

    Uπs = symmetry_adapted_basis(T, tbl, ehom; semisimple = semisimple)

    basis = StarAlgebras.Basis{UInt32}(basis_full)
    invariants = invariant_vectors(tbl, action, basis)

    if try_harder
        Uπs .= simple_projection.(Uπs, Ref(G), Ref(ehom))
    end

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
function Base.eltype(wbdec::WedderburnDecomposition)
    return eltype(eltype(direct_summands(wbdec)))
end

function diagonalize(
    A::AbstractMatrix,
    wbdec::WedderburnDecomposition;
    trace_preserving::Bool = true,
)
    T = promote_type(eltype(A), eltype(wbdec))
    As = [
        (d = size(Uπ, 1); similar(A, T, d, d)) for Uπ in direct_summands(wbdec)
    ]
    return diagonalize!(As, A, wbdec; trace_preserving = trace_preserving)
end

function diagonalize!(
    As::AbstractVector{<:AbstractMatrix},
    A::AbstractMatrix,
    wbdec::WedderburnDecomposition;
    trace_preserving::Bool = true,
)
    dsummands = direct_summands(wbdec)
    @assert axes(As) == axes(dsummands)
    @assert all(size(ds, 1) == size(M, 1) for (ds, M) in zip(dsummands, As))
    @assert all(size(ds, 2) == size(A, 1) for ds in dsummands)
    @assert all(==(size(M)...) for M in As)

    _eps = length(dsummands) * size(A, 1) * eps(eltype(wbdec))

    Threads.@threads for i in eachindex(As)
        ds = dsummands[i]
        U = image_basis(ds)
        # this is the faster version when U are row-based
        # re-test when U move to column-based
        As[i] .= (U * (A * U'))
        if trace_preserving
            As[i] .*= degree(ds)
        end
        if issparse(As[i])
            if eltype(As[i]) <: AbstractFloat
                SparseArrays.droptol!(As[i], _eps)
            else
                SparseArrays.dropzeros!(As[i])
            end
        end
    end
    if eltype(A) <: AbstractFloat
        if trace_preserving && abs(tr(A) - sum(tr, As)) > _eps
            @warn "decomposition did not preserve the trace; check the invariance of A" tr(
                A,
            ) sum(tr, As)
        end
    end
    return As
end

function invariant_vectors(
    tbl::Characters.CharacterTable,
    act::Action,
    basis::StarAlgebras.Basis,
)
    triv_χ =
        Characters.Character{Rational{Int}}(Characters.trivial_character(tbl))
    ehom = ExtensionHomomorphism(act, basis)

    mpr = matrix_projection_irr(ehom, triv_χ)
    mpr, pivots = row_echelon_form!(mpr)
    img = mpr[1:length(pivots), :]

    # TODO:
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

function invariant_vectors(
    tbl::Characters.CharacterTable,
    act::BySignedPermutations,
    basis::StarAlgebras.Basis{T,I},
) where {T,I}
    G = parent(tbl)
    ordG = order(Int, G)
    elts = collect(G)
    orbit = zeros(I, ordG)
    CT = promote_type(coeff_type(act), Rational{Int}) # output coeff type
    coeffs = Vector{CT}(undef, ordG)
    tovisit = trues(length(basis))
    invariant_vs = Vector{SparseVector{CT,Int}}()

    sizehint!(invariant_vs, length(basis) ÷ ordG)

    for (i, b) in enumerate(basis)
        if tovisit[i]
            Threads.@threads for j in eachindex(elts)
                g = elts[j]
                gb, c = SymbolicWedderburn.action(act, g, b)
                orbit[j] = basis[gb]
                coeffs[j] = c
            end
            @view(tovisit[orbit]) .= false
            v = sparsevec(orbit, 1 // ordG .* coeffs, length(basis))
            if (VT = eltype(v)) <: Union{AbstractFloat,Complex}
                droptol!(v, eps(real(VT)) * length(v))
            end
            if !iszero(v)
                push!(invariant_vs, v)
            end
        end
    end
    return invariant_vs
end
