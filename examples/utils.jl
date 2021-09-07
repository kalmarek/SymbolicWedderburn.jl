
using LinearAlgebra
using SparseArrays
using SymbolicWedderburn.StarAlgebras
using GroupsCore

function orbit_constraint!(
    M_orb::Matrix{<:AbstractFloat},
    M::Matrix{<:Integer},
    orb
)
    M_orb .= zero(eltype(M_orb))
    k = 1/length(orb)
    for i in eachindex(M)
        if M[i] ∈ orb
            M_orb[i] += k
        end
    end
    return M_orb
end

# sum(M.==idx for idx in orb)


function orbit_decomposition(
    Γ::Group,
    act::SymbolicWedderburn.ByPermutations,
    basis::StarAlgebras.Basis{T, I},
) where {T, I}

    tovisit = trues(length(basis));
    orbits = Vector{Vector{I}}()

    elts = collect(Γ)
    orbit = zeros(I, length(elts))

    for i in eachindex(basis)
        if tovisit[i]
            g = basis[i]
            Threads.@threads for j in eachindex(elts)
                orbit[j] = basis[SymbolicWedderburn.action(act, elts[j], g)]
            end
            tovisit[orbit] .= false
            push!(orbits, unique(orbit))
        end
    end
    return orbits
end

struct WedderburnDecomposition{B, I, DS}
    basis::B
    orbits::Vector{Vector{I}}
    Uπs::Vector{DS}
    hom::SymbolicWedderburn.InducedActionHomomorphism
end

function WedderburnDecomposition(
    T::Type,
    G::Group,
    action::SymbolicWedderburn.ByPermutations,
    basis_full,
    basis_half,
    S = Rational{Int};
)
    basis = StarAlgebras.Basis{UInt32}(basis_full)
    orbits = orbit_decomposition(G, action, basis)

    tbl = SymbolicWedderburn.CharacterTable(S, G)
    ehom = SymbolicWedderburn.CachedExtensionHomomorphism(G, action, basis_half, precompute=true)

    Uπs = let
        sa_basis = SymbolicWedderburn.symmetry_adapted_basis(T, tbl, ehom; semisimple=false)
        @info round.(_fillfactor.(sa_basis), digits=3)
        sp_basis = sparse.(sa_basis)
        if T <: AbstractFloat
            droptol!.(sp_basis, eps(T)*length(basis_half))
            @info round.(_fillfactor.(sa_basis), digits=3)
        end
        sp_basis
    end

    return WedderburnDecompositionA(basis, orbits, Uπs, ehom)
end

orbits(wbdec::WedderburnDecomposition) = wbdec.orbits
SymbolicWedderburn.basis(wbdec::WedderburnDecomposition) = wbdec.basis
projections(wbdec::WedderburnDecomposition) = wbdec.Uπs

_tmps(wbdec::WedderburnDecomposition) = [similar(U, reverse(size(U))) for U in wbdec.Uπs]

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
    return diagonalize!(Mπs, M, wbdec, tmps)
end

function diagonalize!(
    Mπs,
    M::AbstractMatrix,
    wbdec::WedderburnDecomposition,
    tmps=_tmps(wbdec))

    for (π, Uπ) in pairs(projections(wbdec))
        bUπ = SymbolicWedderburn.basis(Uπ) # Base.Matrix to allow BLAS paths below
        LinearAlgebra.mul!(tmps[π], M, bUπ')
        LinearAlgebra.mul!(Mπs[π], bUπ, tmps[π])
        deg = SymbolicWedderburn.degree(Uπ)
        Mπs[π] .*= deg
    end
    return Mπs
end

