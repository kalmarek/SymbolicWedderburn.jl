"""
    matrix_projection(χ::AbstractClassFunction)
Compute matrix projection associated to the abstract class function `χ`.

The dimension `d` of the projection is derived from `conjugacy_classes(χ)`.
E.g. if `conjugacy_classes(χ)` consist of permutations of degree `d` (i.e.
acting naturally on the set `1:d`) the result will be a matrix of size `(d,d)`.

Returned tuple consist of the matrix realization of the projection and a coefficient (its weight).
"""
function matrix_projection(χ::Character{T}) where {T}
    U = matrix_projection(values(χ), conjugacy_classes(χ))
    deg = degree(χ)
    ordG = sum(length, conjugacy_classes(χ))

    return U, T(deg) / ordG
end

function matrix_projection(χ::AbstractClassFunction)
    deg = degree(χ)
    if iszero(deg) # short circuting the trivial case
        return zeros(eltype(χ), 0, degree(first(first(conjugacy_classes(χ))))), deg
    end
    ordG = sum(length, conjugacy_classes(χ))
    U = matrix_projection(values(χ), conjugacy_classes(χ))

    return U, deg / ordG
end

"""
    matrix_projection(vals, ccls)
Return matrix projection defined by an abstract class function with values `vals`
attained on (the corresponding) conjugacy classes `ccls`.

The dimension of the projection is equal to the degree of the permutations in `ccls`.
"""
function matrix_projection(
    vals::AbstractVector{T},
    ccls::AbstractVector{<:AbstractOrbit{<:AbstractPerm}},
    dim = degree(first(first(ccls))),
) where {T}

    result = zeros(T, dim, dim)

    for (val, cc) in zip(vals, ccls)
        for g in cc
            for i in 1:dim
                result[i, i^g] += val
            end
        end
    end

    return result
end

function matrix_projection(
    vals::AbstractVector{T},
    ccls::AbstractVector{<:AbstractOrbit{<:AbstractMatrix}},
    dim = size(first(first(ccls)), 1),
) where {T}

    return sum(val .* sum(g -> inv(Matrix(g)), cc) for (val, cc) in zip(vals, ccls))
end

"""
    action_character(conjugacy_cls)
Return the character of the representation given by the elements in the
conjugacy classes `conjugacy_cls`.

This corresponds to the classical definition of characters as a traces of the
representation matrices.
"""
function action_character(
    conjugacy_cls::AbstractVector{CCl},
    inv_of = _inv_of(conjugacy_cls),
) where {CCl<:AbstractOrbit{<:AbstractPerm}}
    vals = Int[PermutationGroups.nfixedpoints(first(cc)) for cc in conjugacy_cls]
    # in general:
    # vals = [tr(matrix_representative(first(cc))) for cc in conjugacy_cls]
    return Character(vals, inv_of, conjugacy_cls)
end

function action_character(
    conjugacy_cls::AbstractVector{CCl},
    inv_of = _inv_of(conjugacy_cls),
) where {CCl<:AbstractOrbit{<:AbstractMatrix}}
    vals = [tr(first(cc)) for cc in conjugacy_cls]
    return Character(vals, inv_of, conjugacy_cls)
end

"""
    affordable_real!(chars::AbstractVector{<:AbstractClassFunction})
Return _real_ characters formed from `chars` by replacing `χ` with `2re(χ)` when necessary.
"""
function affordable_real!(chars::AbstractVector{T}) where {T<:AbstractClassFunction}
    all(all(isreal, values(χ)) for χ in chars) && return chars
    pmap = PowerMap(conjugacy_classes(first(chars)))
    res = affordable_real!.(chars, Ref(pmap))
    return unique!(res)
end
