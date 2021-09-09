
using LinearAlgebra
using SparseArrays
using SymbolicWedderburn.StarAlgebras
using GroupsCore

function invariant_constraint!(
    M_orb::AbstractMatrix{<:AbstractFloat},
    M::Matrix{<:Integer},
    invariant_vec::SparseVector,
)
    M_orb .= zero(eltype(M_orb))
    for i in eachindex(M)
        if M[i] ∈ SparseArrays.nonzeroinds(invariant_vec)
            M_orb[i] += invariant_vec[M[i]]
        end
    end
    return M_orb
end

function invariant_vectors(
    G::Group,
    act::SymbolicWedderburn.ByPermutations,
    basis::StarAlgebras.Basis{T, I},
) where {T, I}

    tovisit = trues(length(basis));
    invariant_vs = Vector{SparseVector{Rational{Int}}}()

    ordG = order(Int, G)
    elts = collect(G)
    orbit = zeros(I, ordG)

    for i in eachindex(basis)
        if tovisit[i]
            g = basis[i]
            Threads.@threads for j in eachindex(elts)
                orbit[j] = basis[SymbolicWedderburn.action(act, elts[j], g)]
            end
            tovisit[orbit] .= false
            push!(invariant_vs, sparsevec(orbit, fill(1//ordG, ordG), length(basis)))
        end
    end
    return invariant_vs
end

_orth_proj(A::AbstractMatrix, v, QR=qr(A)) = A * (QR\v)

function invariant_vectors(
    G::Group,
    act::SymbolicWedderburn.ByLinearTransformation,
    basis,
    atol=1e-12
    )

    hom = SymbolicWedderburn.CachedExtensionHomomorphism(G, act, basis, precompute = true)

    dim = 0
    T = SymbolicWedderburn.coeff_type(SymbolicWedderburn.action(hom))
    @assert T==Float64
    subspaces = SparseMatrixCSC{T, Int64}[]
    qrs = LinearAlgebra.QRCompactWY{Float64, Matrix{Float64}}[]

    inv_vectors = SparseVector{T, Int64}[]
    inv_v = Vector{T}(undef, length(basis))

    for b in basis
        b_orth = sparsevec(SymbolicWedderburn.decompose(b, hom)..., length(basis))
        if dim > 0
            b_orth -= sum(_orth_proj(V, b_orth, QR) for (V, QR) in zip(subspaces, qrs))
            @debug "orthogonalizing $b:" norm(b_orth)
            norm(b_orth, 2) < atol && continue
        end

        # we found a new subspace!
        b_orth ./= norm(b_orth, 2)
        orbit_subspace = Matrix{T}(undef, length(basis), 1)
        orbit_subspace[:, 1] .= b_orth
        inv_v .= zero(T)

        # for the orbit of b_orth find
        # * `inv_v`, the invariant (average) vector, and
        # * `tmp_subspace`, the subspace spanned by the orbit

        for g in G
            k = SymbolicWedderburn.induce(hom, g) * b_orth
            @assert isapprox(norm(k, 2), 1.0, atol=atol)
            inv_v .+= k

            residual = k - _orth_proj(orbit_subspace, k)
            (rnorm = norm(residual, 2)) < atol && continue
            orbit_subspace = hcat(orbit_subspace, Vector(residual./rnorm))
        end

        push!(inv_vectors, inv_v/order(Int, G)) # a copy if inv_v
        push!(subspaces, orbit_subspace) # orbit_subspace is freshly allocated every nonzero b_orth
        dim += size(orbit_subspace, 2)
        dim == length(basis) && break
        dim > length(basis) && throw("This should never occur")
        # moving qr here to save one qr decomposition
        push!(qrs, qr(orbit_subspace))
    end
    @assert dim == length(basis) "Encountered incomplete decomposition into orbits!"
    return inv_vectors
end


struct WedderburnDecomposition{B, O, DS<:SymbolicWedderburn.DirectSummand}
    basis::B
    invariants::O
    Uπs::Vector{DS}
    hom::SymbolicWedderburn.InducedActionHomomorphism
end

function WedderburnDecomposition(
    T::Type,
    G::Group,
    action::SymbolicWedderburn.Action,
    basis_full,
    basis_half,
    S = Rational{Int};
)
    basis = StarAlgebras.Basis{UInt32}(basis_full)
    invariants = invariant_vectors(G, action, basis)

    tbl = SymbolicWedderburn.CharacterTable(S, G)
    ehom = SymbolicWedderburn.CachedExtensionHomomorphism(G, action, basis_half, precompute=true)

    Uπs = let
        sa_basis = SymbolicWedderburn.symmetry_adapted_basis(T, tbl, ehom; semisimple=false)
        @debug "fill factor of sa_basis:" round.(_fillfactor.(sa_basis), digits=3)
        sp_basis = sparse.(sa_basis)
        if T <: AbstractFloat
            droptol!.(sp_basis, eps(T)*length(basis_half))
            @debug "fill factor after sparsification" round.(_fillfactor.(sa_basis), digits=3)
        end
        sp_basis
    end

    return WedderburnDecomposition(basis, invariants, Uπs, ehom)
end

invariant_vectors(wbdec::WedderburnDecomposition) = wbdec.invariants
SymbolicWedderburn.basis(wbdec::WedderburnDecomposition) = wbdec.basis
projections(wbdec::WedderburnDecomposition) = wbdec.Uπs

_tmps(wbdec::WedderburnDecomposition) =
    [zeros(eltype(U), reverse(size(U))) for U in wbdec.Uπs]

_fillfactor(M::AbstractMatrix) = count(!iszero, M)/length(M)
_fillfactor(M::AbstractSparseMatrix) = nnz(M)/length(M)

function diagonalize(
    M::AbstractMatrix,
    wbdec::WedderburnDecomposition,
    tmps=_tmps(wbdec)
)
    # return [SymbolicWedderburn.degree(U).*(U*(M*U')) for U in wbdec.Uπs]

    T = eltype(eltype(projections(wbdec)))
    Mπs = [(d = size(Uπ, 1); zeros(T, d,d)) for Uπ in projections(wbdec)]
    return diagonalize!(Mπs, M, projections(wbdec), tmps)
end

function diagonalize!(
    Mπs,
    M::AbstractMatrix,
    wbdec::WedderburnDecomposition,
    tmps=_tmps(wbdec))
    return diagonalize!(Mπs, M, projections(wbdec), tmps)
end

function diagonalize!(
    Mπs::AbstractVector{<:AbstractMatrix},
    M::AbstractMatrix,
    Uπs::AbstractVector{<:SymbolicWedderburn.DirectSummand},
    tmps
)
    for (π, Uπ) in enumerate(Uπs)
        bUπ = SymbolicWedderburn.basis(Uπ) # Base.Matrix to allow BLAS paths below
        LinearAlgebra.mul!(tmps[π], M, bUπ')
        LinearAlgebra.mul!(Mπs[π], bUπ, tmps[π])
        deg = SymbolicWedderburn.degree(Uπ)
        Mπs[π] .*= deg
    end
    return Mπs
end

_diagonalize(M, Uπs) =[SymbolicWedderburn.degree(U).*(U*(M*U')) for U in Uπs]
