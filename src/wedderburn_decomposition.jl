struct WedderburnDecomposition{B,iV,DS<:DirectSummand,Hom}
    basis::B
    invariants::iV
    Uπs::Vector{DS}
    hom::Hom
end

"""
    WedderburnDecomposition(T::Type, G::Group, action::Action, basis_full, basis_half, S = Rational{Int}; semisimple = false)
Compute a Wedderburn decomposition of for group *-algebras arising from (non)commutative sum-of-squares optimization.

Return a WedderburnDecomposition consisting of
- basis: the original basis_full
- invariants: a basis for the invariant subspace <basis_full>^G
- Uπs: a Vector (indexed by irrep type i of G) of `DirectSummand`s of the G-modules `V=span(basis_half)`
- hom: stores elements of G as endomorphism of V (a sparse matrix if action<:ByLinearTransformation, or a permutation if action<:ByPermutations)


Remaining arguments (T, S, semisimple) control the behavior of the algorithm and its output
- basis: the original basis_full
- invariants: a basis of End(V)^G as a linear space
- Uπs: a vector (indexed by irrep type i of G) of `DirectSummand`s of End_G(V)
- hom: extension homomorphism ???
Each summand `Uπs[i]` contains the following fields:
- basis: a **dense** mᵢnᵢ×length(basis_half) [depends on semisimple: if true, correct] matrix giving the
orthogonal projection to this isotypic summand. That is X ↦ Uπ.basis X Uπ.basis^† projects
X ∈ End(V) to the summand id_Vᵢ⊗M_mᵢ(k) [in V, projects to mᵢVᵢ and adjoint to id_Vᵢ⊗M_mᵢ(k) as a direct sum]
- degree: dimension of an irrep of this type
- multiplicity: multiplicity of this irrep in this isotypic component
- simple: true if `Wᵢ≅mᵢVᵢ` (given by basis) is an isomorphism

The

See also: [symmetry_adapted_basis](@ref)
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
    basis = StarAlgebras.Basis{UInt32}(basis_full)
    invariants = invariant_vectors(G, action, basis)

    tbl = CharacterTable(S, G)
    ehom = CachedExtensionHomomorphism(G, action, basis_half, precompute = true)

    Uπs = let
        sa_basis = symmetry_adapted_basis(T, tbl, ehom; semisimple = semisimple)
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
    G::Group,
    act::ByPermutations,
    basis::StarAlgebras.Basis{T,I},
) where {T,I}
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

function _orth_proj(A::AbstractMatrix, v, QR = qr(A))
    # return A * (QR \ v)
    y = QR.Q' * v
    y[size(A, 2)+1:end] .= 0 # project to y = Q̂'x
    return QR.Q * y
end

function invariant_vectors(
    G::Group,
    act::ByLinearTransformation,
    basis,
    atol = 1e-12,
)
    hom = CachedExtensionHomomorphism(G, act, basis, precompute = true)

    dim = 0
    T = coeff_type(action(hom))
    @assert T == Float64
    subspaces = SparseMatrixCSC{T,Int64}[]
    qrs = LinearAlgebra.QRCompactWY{Float64,Matrix{Float64}}[]

    inv_vectors = SparseVector{T,Int64}[]
    inv_v = Vector{T}(undef, length(basis))

    for b in basis
        b_orth = sparsevec(decompose(b, hom)..., length(basis))
        if dim > 0
            b_orth -= sum(
                _orth_proj(V, b_orth, QR) for (V, QR) in zip(subspaces, qrs)
            )
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
            k = induce(hom, g) * b_orth
            @assert isapprox(norm(k, 2), 1.0, atol = atol)
            inv_v .+= k

            residual = k - _orth_proj(orbit_subspace, k)
            (rnorm = norm(residual, 2)) < atol && continue
            orbit_subspace = hcat(orbit_subspace, Vector(residual ./ rnorm))
        end

        push!(inv_vectors, inv_v / order(Int, G)) # a copy if inv_v
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
