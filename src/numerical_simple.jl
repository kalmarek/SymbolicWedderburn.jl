# Numerical fallback used by `_symmetry_adapted_basis` when the symbolic
# `minimal_projection_system` cannot reduce an isotypical component to a
# simple summand (i.e. the returned projection has rank > 1). For
# `T<:BlasFloat`, we compute the matrix representation of the group on the
# isotypical subspace via `induce(hom, g)`, simultaneously block-diagonalize
# via Schur decomposition, and extract one of the d identical simple blocks.

is_orthogonal(::Any) = false
function is_orthogonal(S::AbstractMatrix, tol = 1e-6)
    n = LinearAlgebra.checksquare(S)
    for i in 1:n
        for j in 1:n
            if i != j
                s = LinearAlgebra.dot(view(S, :, i), view(S, :, j))
                if abs(s) > tol
                    return false
                end
            end
        end
    end
    return true
end

_isapproxless(a::Real, b::Real) = a < b
function _isapproxless(a::Complex, b::Complex)
    if real(a) ≈ real(b)
        return isless(imag(a), imag(b))
    else
        return isless(real(a), real(b))
    end
end

"""
    _reorder!(F::LinearAlgebra.Schur{T}) where {T}

Reorder a Schur decomposition so that eigenvalues appear in increasing order.
For `T<:Real`, complex pairs are kept together.
"""
function _reorder!(F::LinearAlgebra.Schur{T}) where {T}
    n = length(F.values)
    sorted = false
    while !sorted
        prev_i = nothing
        sorted = true
        i = 1
        while i <= n
            S = F.Schur
            if (T <: Real) && i < n && !iszero(S[i+1, i])
                next_i = i + 2
            else
                next_i = i + 1
            end
            if !isnothing(prev_i) && _isapproxless(S[i, i], S[prev_i, prev_i])
                select = trues(n)
                select[prev_i:(i-1)] .= false
                select[next_i:end] .= false
                LinearAlgebra.ordschur!(F, select)
                sorted = false
            end
            prev_i = i
            i = next_i
        end
    end
end

# Diagonal `d` with `d[i]*conj(d[i]) = 1` (so `±1` in the real case, unit
# complex in the complex case). Chooses signs/phases so that `Ai` and `Bi`
# agree on the strict upper-triangular part.
function _sign_diag(
    As::Vector{<:AbstractMatrix{T}},
    Bs::Vector{<:AbstractMatrix{T}};
    tol = Base.rtoldefault(real(T)),
) where {T}
    n = LinearAlgebra.checksquare(As[1])
    d = ones(T, n)
    for j in 2:n
        if T <: Real
            minus = zero(real(T))
            not_minus = zero(real(T))
            for i in 1:(j-1)
                for (Ai, Bi) in zip(As, Bs)
                    a = Ai[i, j]
                    b = Bi[i, j]
                    minus = max(minus, abs(a + b))
                    not_minus = max(not_minus, abs(a - b))
                end
            end
            if minus < not_minus
                d[j] = -one(T)
                for B in Bs
                    B[:, j] = -B[:, j]
                    B[j, :] = -B[j, :]
                end
            end
        else
            k = argmax(eachindex(Bs)) do k
                return maximum(abs.(Bs[k][1:(j-1), j]))
            end
            i = argmax(abs.(Bs[k][1:(j-1), j]))
            if abs(Bs[k][i, j]) <= tol
                continue
            end
            rot = As[k][i, j] / Bs[k][i, j]
            rot /= abs(rot)
            d[j] = rot
            for B in Bs
                B[:, j] *= rot
                B[j, :] *= conj(rot)
            end
        end
    end
    return d
end

# Snap to integers if every entry is within rounding of one.
function _try_integer!(A::Matrix)
    if all(a -> round(a) ≈ a, A)
        return round.(A)
    else
        return A
    end
end

"""
    _rotate_complex(A, B; tol)

Given (quasi) upper-triangular `A` and `B` whose real eigenvalues already align
and whose complex pairs may differ by sign/order, return a signed permutation
`P` such that `P'*A*P` and `B` match on the strict-lower-triangular part.
"""
function _rotate_complex(
    A::AbstractMatrix{T},
    B::AbstractMatrix{T};
    tol = Base.rtoldefault(real(T)),
) where {T}
    n = LinearAlgebra.checksquare(A)
    I = collect(1:n)
    J = copy(I)
    V = ones(T, n)
    pair = false
    for i in 1:n
        if pair || i == n
            continue
        end
        pair = abs(A[i+1, i]) > tol
        if pair
            a = (A[i+1, i], A[i, i+1])
            b = (B[i+1, i], B[i, i+1])
            c = a[2:-1:1]
            if LinearAlgebra.norm(abs.(a) .- abs.(b)) >
               LinearAlgebra.norm(abs.(c) .- abs.(b))
                a = c
                J[i] = i + 1
                J[i+1] = i
            end
            c = (-).(a)
            if LinearAlgebra.norm(a .- b) > LinearAlgebra.norm(c .- b)
                V[i+1] = -V[i]
            end
        end
    end
    return SparseArrays.sparse(I, J, V, n, n)
end

"""
    orthogonal_transformation_to(A, B, Ais, Bis)

Return an orthogonal `U` with `A ≈ U'*B*U` and `Ai[k] ≈ U'*Bi[k]*U` for all
pairs in `zip(Ais, Bis)`. Used to align two equivalent diagonal blocks.
"""
function orthogonal_transformation_to(A, B, Ais, Bis)
    As = LinearAlgebra.schur(A)
    _reorder!(As)
    T_A = As.Schur
    Z_A = As.vectors
    Bs = LinearAlgebra.schur(B)
    _reorder!(Bs)
    T_B = Bs.Schur
    Z_B = Bs.vectors
    P = _rotate_complex(T_A, T_B)
    T_A = P' * T_A * P
    Ais = [P' * Z_A' * A * Z_A * P for A in Ais]
    push!(Ais, T_A)
    Bis = [Z_B' * B * Z_B for B in Bis]
    push!(Bis, T_B)
    d = _sign_diag(Ais, Bis)
    D = LinearAlgebra.Diagonal(d)
    return _try_integer!(Z_B * D * P' * Z_A')
end

# Minimal union-find used to detect block structure in the Schur basis.
mutable struct _UnionFind
    parents::Vector{Int}
    ranks::Vector{Int}
    _UnionFind(n::Integer) = new(collect(1:n), zeros(Int, n))
end

function _find_root!(uf::_UnionFind, i::Int)
    p = uf.parents[i]
    if p == i
        return i
    end
    r = _find_root!(uf, p)
    uf.parents[i] = r
    return r
end

function _union!(uf::_UnionFind, i::Int, j::Int)
    ri = _find_root!(uf, i)
    rj = _find_root!(uf, j)
    ri == rj && return ri
    if uf.ranks[ri] < uf.ranks[rj]
        uf.parents[ri] = rj
        return rj
    elseif uf.ranks[ri] > uf.ranks[rj]
        uf.parents[rj] = ri
        return ri
    else
        uf.parents[rj] = ri
        uf.ranks[ri] += 1
        return ri
    end
end

function _merge_sparsity!(uf::_UnionFind, A, tol = 1e-8)
    for I in CartesianIndices(A)
        i, j = I.I
        if abs(A[I]) > tol
            _union!(uf, i, j)
        end
    end
end

function _isblockdim(A, d, tol = 1e-8)
    for I in CartesianIndices(A)
        i, j = I.I
        if abs(i - j) >= d && abs(A[I]) > tol
            return false
        end
    end
    return true
end

function _is_ordered_blockdim(A, d, tol = 1e-8)
    B = A[1:d, 1:d]
    return _isblockdim(A, d, tol) && all(d:d:(size(A, 1)-d)) do offset
        I = offset .+ (1:d)
        return isapprox(B, A[I, I], rtol = tol)
    end
end

function _ordered_block_check(U, As, d)
    iU = U'
    if !(iU ≈ inv(U))
        return false
    end
    return all(As) do A
        return _is_ordered_blockdim(iU * A * U, d)
    end
end

# Block-diagonalize: find an orthogonal `U` such that for every `A ∈ As`,
# `U' A U` is block-diagonal whose blocks are multiples of `d` in size.
function block_diagonalize(As, d)
    for A in As
        Z = LinearAlgebra.schur(A).vectors
        iZ = Z'
        @assert iZ ≈ inv(Z)
        n = LinearAlgebra.checksquare(A)
        uf = _UnionFind(n)
        Bs = [iZ * A * Z for A in As]
        for B in Bs
            _merge_sparsity!(uf, B)
        end
        blocks = Dict{Int,Vector{Int}}()
        for i in 1:n
            r = _find_root!(uf, i)
            blocks[r] = push!(get(blocks, r, Int[]), i)
        end
        if length(blocks) > 1
            U = similar(Z)
            offset = 0
            broke = false
            for v in values(blocks)
                @assert iszero(length(v) % d)
                V = Z[:, v]
                if length(v) != d
                    Cs = [B[v, v] for B in Bs]
                    _V = block_diagonalize(Cs, d)
                    if isnothing(_V)
                        broke = true
                        break
                    end
                    V *= _V
                end
                U[:, offset .+ eachindex(v)] = V
                offset += length(v)
            end
            broke && continue
            if offset == n
                return U
            end
        end
    end
    return nothing
end

"""
    ordered_block_diagonalize(As, d)

Return an orthogonal `U` with `U'*A*U` block-diagonal in `d×d` blocks for every
`A ∈ As`, with **all** blocks equal to a common reference. Returns `nothing` if
no such `U` exists.
"""
function ordered_block_diagonalize(As, d)
    U = block_diagonalize(As, d)
    isnothing(U) && return nothing
    iU = U'
    @assert iU ≈ inv(U)
    Bs = [iU * A * U for A in As]
    @assert all(B -> _isblockdim(B, d), Bs)
    refs = [B[1:d, 1:d] for B in Bs]
    for offset in d:d:(size(U, 1)-d)
        I = offset .+ (1:d)
        Cs = [B[I, I] for B in Bs]
        # Random combination trick (cf. [CGT97]) to align blocks with probability 1.
        λ = rand(length(Bs))
        R = sum(λ .* refs)
        C = sum(λ .* Cs)
        V = orthogonal_transformation_to(R, C, refs, Cs)
        @assert R ≈ V' * C * V
        for i in eachindex(refs)
            @assert refs[i] ≈ V' * Cs[i] * V
        end
        U[:, I] = U[:, I] * V
    end
    @assert _ordered_block_check(U, As, d)
    return U
end

# Build matrix representations of the group generators restricted to the
# isotypical subspace spanned by the rows of `R` (size `(m*d, N)`), using the
# induced action homomorphism.
#
# For orthonormal-rowed `R` and orthogonal action matrices `M_g = induce(hom, g)`,
# the matrix `S_g = R * M_g * R'` represents `g` on the basis given by the
# rows of `R`. Block structure is invariant under transpose, so even if a
# convention mismatch flips `S_g ↔ S_g'`, `ordered_block_diagonalize` still
# returns the correct `U`.
function _isotypical_matrix_reps(
    R::AbstractMatrix{T},
    hom::InducedActionHomomorphism,
    G,
) where {T}
    n = size(R, 2)
    return map(GroupsCore.gens(G)) do g
        M = _dense_action_matrix(T, induce(hom, g), n)
        return Matrix{T}(R * M * transpose(R))
    end
end

# `induce(hom, g)` returns either a permutation or a sparse matrix; densify.
function _dense_action_matrix(::Type{T}, M::AbstractMatrix, _) where {T}
    return convert(Matrix{T}, M)
end
function _dense_action_matrix(
    ::Type{T},
    p::PG.AbstractPermutation,
    n::Integer,
) where {T}
    A = zeros(T, n, n)
    for i in 1:n
        A[i, i^p] = one(T)
    end
    return A
end

"""
    numerical_simplify(ds::DirectSummand, hom::InducedActionHomomorphism, G)

Numerically reduce a non-simple `DirectSummand` to a simple one of shape
`m × N` by simultaneously block-diagonalizing the matrix representations of
the generators of `G` on the isotypical subspace spanned by `image_basis(ds)`.

If `ds` is already simple this is a no-op. Throws if the induced action on
the chosen image basis is not orthogonal (in which case the symbolic
decomposition should be preferred).

Requires `eltype(image_basis(ds)) <: BlasFloat` for the Schur decomposition.
"""
function numerical_simplify(
    ds::DirectSummand,
    hom::InducedActionHomomorphism,
    G,
)
    issimple(ds) && return ds
    R = image_basis(ds)
    T = eltype(R)
    T <: LinearAlgebra.BlasFloat || throw(ArgumentError(
        "numerical_simplify requires `eltype(image_basis(ds)) <: BlasFloat`, got $T",
    ))
    m = multiplicity(ds)
    d = AP.degree(character(ds))
    # `R` should span the full m*d-dim isotypical subspace. The caller is
    # responsible for re-projecting onto rank m*d before calling this.
    @assert size(R, 1) == m * d "`numerical_simplify` expects R with m*d rows; got $(size(R, 1)) for m=$m d=$d"

    F = convert(Matrix{T}, R)
    Ss = _isotypical_matrix_reps(Matrix{T}(R), hom, G)
    if !all(is_orthogonal, Ss)
        error("The matrix representation induced from the action on the polynomial basis is not orthogonal.")
    end
    U = ordered_block_diagonalize(Ss, d)
    if isnothing(U)
        error("Could not simultaneously block-diagonalize into $m identical $(d)x$(d) blocks")
    end
    # Pick one of the d simple sub-blocks. For SDP applications, all d blocks
    # are equivalent (give the same Gram matrix `C`); see `diagonalize!` which
    # accounts for the d-fold replication via `trace_preserving`.
    cols = 1:d:(1 + d*(m-1))
    simple = transpose(U[:, cols]) * F
    return DirectSummand(simple, m, character(ds))
end
