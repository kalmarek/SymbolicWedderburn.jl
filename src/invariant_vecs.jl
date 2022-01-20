struct InvariantVectors{T,A}
    rows::A
    size::NTuple{2,Int}
end
InvariantVectors(M::AbstractMatrix) = InvariantVectors{eltype(M)}(M)
function InvariantVectors{T}(M::AbstractMatrix) where {T}
    m = convert.(T, M)
    return InvariantVectors{T,typeof(m)}(m, size(m))
end

InvariantVectors(rows::AbstractVector{<:AbstractVector{T}}) where {T} =
    InvariantVectors{T}(rows)

function InvariantVectors{T}(
    rows::AbstractVector{<:AbstractVector},
    n = length(first(rows))
) where {T}
    @assert all(==(n) ∘ length, rows)

    a = broadcast.(x -> convert(T, x), rows)
    return InvariantVectors{T,typeof(a)}(a, (length(a), n))
end

Base.length(iv::InvariantVectors) = first(iv.size)
Base.eltype(iv::InvariantVectors{T, <:AbstractVector{S}}) where {T, S} = S
Base.eltype(iv::InvariantVectors{T, M}) where {T, M<:AbstractMatrix} =
    SubArray{T, 1, M, Tuple{Int64, Base.Slice{Base.OneTo{Int64}}}, true}

Base.iterate(iv::InvariantVectors) = iterate(iv.rows)
Base.iterate(iv::InvariantVectors, state) = iterate(iv.rows, state)

function Base.iterate(iv::InvariantVectors{T, M}) where {T, M<:AbstractMatrix}
    isempty(iv.rows) && return nothing
    itr = eachrow(iv.rows)
    row, st = iterate(itr)
    return row, (itr, st)
end

function Base.iterate(iv::InvariantVectors{T, M}, state) where {T, M<:AbstractMatrix}
    itr, st = state
    (tmp = iterate(itr, st)) === nothing && return nothing
    row, st = tmp
    return row, (itr, st)
end

function invariant_vectors(
    S::Type,
    tbl::Characters.CharacterTable,
    act::Action,
    basis::StarAlgebras.Basis,
)
    triv_χ =
        Characters.Character{Rational{Int}}(Characters.trivial_character(tbl))
    ehom =
        CachedExtensionHomomorphism(parent(tbl), act, basis, precompute = true)
    # ehom = ExtensionHomomorphism(act, basis)

    mpr = matrix_projection_irr(ehom, triv_χ)
    mpr, pivots = row_echelon_form!(mpr)
    img = mpr[1:length(pivots), :]

    return InvariantVectors{S}(img)
end

function invariant_vectors(
    S::Type,
    tbl::Characters.CharacterTable,
    act::ByPermutations,
    basis::StarAlgebras.Basis{T,I},
) where {T,I}
    G = parent(tbl)
    tovisit = trues(length(basis))
    invariant_vs = Vector{SparseVector{S}}()

    elts = collect(G)
    ordG = length(elts)
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
    return InvariantVectors{S}(invariant_vs)
end
