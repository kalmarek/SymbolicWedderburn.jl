"""
    AbstractClassFunction{T}
Abstract type representing functions constant on conjugacy classes of a group.

The following functionality is required for an `AbstractClassFunction`:
 * `parent(χ)`: the underlying group
 * `conjugacy_classes(χ)`: iterator over conjugacy classes of `χ`.
 * `values(χ)`: a _typed_ iterator over values
 * `getindex(χ, i::Integer)`: the value on `i`-th conjugacy class.
 Note: Indexing with negative integers should return values on the class which
 contains inverses of the `i`-th class.
 * `χ(g::GroupElement)`: value of `χ` on a group element (only for julia < 1.3)
"""
abstract type AbstractClassFunction{T} end # <:AbstractVector{T}?

Base.eltype(::AbstractClassFunction{T}) where T = T

function LinearAlgebra.dot(χ::AbstractClassFunction, ψ::AbstractClassFunction)
    @assert conjugacy_classes(χ) === conjugacy_classes(ψ)
    val = sum(
        length(cc) * χ[i] * ψ[-i] for (i, cc) in enumerate(conjugacy_classes(χ))
    )
    orderG = sum(length, conjugacy_classes(χ))
    val = _div(val, orderG)
    return val
end

_div(val, orderG) = div(val, orderG)
_div(val::AbstractFloat, orderG) = val / orderG
_div(val::ComplexF64, orderG) = val / orderG

## Arbitrary ClassFunctions without decomposition into irreps
struct ClassFunction{T,CCl} <: AbstractClassFunction{T}
    vals::Vector{T}
    conjugacy_classes::CCl
    inv_of::Vector{Int}
end

ClassFunction(vals, cclasses) = ClassFunction(vals, cclasses, _inv_of(cclasses))

## AbstractClassFunction api
Base.parent(χ::ClassFunction) = parent(first(first(conjugacy_classes(χ))))
conjugacy_classes(χ::ClassFunction) = χ.conjugacy_classes
Base.values(χ::ClassFunction) = χ.vals

Base.@propagate_inbounds function Base.getindex(χ::ClassFunction, i::Integer)
    i = i < 0 ? χ.inv_of[-i] : i
    @boundscheck 1 ≤ i ≤ length(conjugacy_classes(χ))
    return values(χ)[i]
end

# TODO: move to AbstractClassFunction when we drop julia-1.0
# adding call methods to abstract types was not supported before julia-1.3
function (χ::ClassFunction)(g::GroupElement)
    for (i, cc) in enumerate(conjugacy_classes(χ))
        g ∈ cc && return χ[i]
    end
    throw(
        DomainError(
            g,
            "element does not belong to conjugacy classes of $χ",
        ),
    )
end

"""
    Character <: AbstractClassFunction
Struct representing (possibly virtual) character of a group.

Characters are backed by `table(χ)::CharacterTable` which actually stores the
character values. The constituents (decomposition into the irreducible summands)
of a given character can be obtained by calling `constituents(χ)` which returns a
vector of coefficients of `χ` in the basis of `irreducible_characters(table(χ))`.

It is assumed that equal class functions on the same group will have **identical**
(ie. `===`) character tables.
"""
struct Character{T, S, ChT<:CharacterTable} <: AbstractClassFunction{T}
    table::ChT
    constituents::Vector{S}
end

function Character(
    chtbl::CharacterTable{Gr, T},
    constituents::AbstractVector{S}
) where {Gr, T, S}
    R = Base._return_type(*, Tuple{S,T})
    return Character{R, S, typeof(chtbl)}(chtbl, constituents)
end

Character(chtbl::CharacterTable, i::Integer) = Character{eltype(chtbl)}(chtbl, i)

function Character{T}(chtbl::CharacterTable, i::Integer) where T
    v = zeros(Int, nconjugacy_classes(chtbl))
    v[i] = 1
    return Character{T, Int, typeof(chtbl)}(chtbl, v)
end

function Character{T}(χ::Character) where T
    S = eltype(constituents(χ))
    ChT = typeof(table(χ))
    return Character{T, S, ChT}(table(χ), constituents(χ))
end

## Accessors
table(χ::Character) = χ.table
constituents(χ::Character) = χ.constituents

## AbstractClassFunction api
Base.parent(χ::Character) = parent(table(χ))
conjugacy_classes(χ::Character) = conjugacy_classes(table(χ))
Base.values(χ::Character) = (χ[i] for i in 1:nconjugacy_classes(table(χ)))

Base.@propagate_inbounds function Base.getindex(χ::Character{T}, i::Integer) where T
    i = i < 0 ? table(χ).inv_of[abs(i)] : i
    @boundscheck 1 ≤ i ≤ nconjugacy_classes(table(χ))
    return T(
        sum(c * table(χ)[idx, i] for (idx, c) in enumerate(constituents(χ))
        if !iszero(c))
    )
end

# TODO: move to AbstractClassFunction when we drop julia-1.0
function (χ::Character)(g::GroupElement)
    for (i, cc) in enumerate(conjugacy_classes(χ))
        g ∈ cc && return χ[i]
    end
    throw(
        DomainError(
            g,
            "element does not belong to conjugacy classes of $χ",
        ),
    )
end

## Basic functionality

Base.:(==)(χ::Character, ψ::Character) =
    table(χ) === table(ψ) && constituents(χ) == constituents(ψ)
Base.hash(χ::Character, h::UInt) = hash(table(χ), hash(constituents(χ), h))

Base.deepcopy_internal(χ::Character, ::IdDict) =
    Character(table(χ), copy(constituents(χ)))

## Character arithmetic

for f in (:+, :-)
    @eval begin
        function Base.$f(χ::Character, ψ::Character)
            @assert table(χ) === table(ψ)
            return Character(table(χ), $f(constituents(χ), constituents(ψ)))
        end
    end
end

Base.:*(χ::Character, c::Number) = Character(table(χ), c .* constituents(χ))
Base.:*(c::Number, χ::Character) = χ * c
Base.:/(χ::Character, c::Number) = Character(table(χ), constituents(χ) ./ c)

## Group-theoretic functions:

PermutationGroups.degree(χ::Character) = Int(χ(one(parent(χ))))
PermutationGroups.degree(
    χ::Character{T,CCl},
) where {T,CCl <: AbstractOrbit{<:AbstractMatrix}} = Int(χ[1])

Base.conj(χ::Character) = Character(table(χ), constituents(χ)[table(χ)._inv_of])

isvirtual(χ::Character) =
    any(<(0), constituents(χ)) || any(!isinteger, constituents(χ))

function isirreducible(χ::Character)
    count(isone, constituents(χ)) == 1 || return false
    count(iszero, constituents(χ)) == length(constituents(χ)) - 1 || return false
    return true
end

"""
    affordable_real!(χ::Character)
Return either `χ` or `2re(χ)` depending whether `χ` is afforded by a real representation, modifying `χ` in place.
"""
function affordable_real!(χ::Character)
    ι = frobenius_schur_indicator(χ)
    if ι <= 0 # i.e. χ is complex or quaternionic
        χ.constituents += constituents(conj(χ))
        return χ
    else
        return χ
    end
end

"""
    frobenius_schur_indicator(χ::AbstractClassFunction[, pmap::PowerMap])
Return Frobenius-Schur indicator of `χ`, i.e. `Σχ(g²)` where sum is taken over
the whole group.

If χ is an irreducible `Character`, Frobenius-Schur indicator takes values in
`{1, 0, -1}` which correspond to the following situations:
 1. `χ` is real-valued and is afforded by an irreducible real representation,
 2. `χ` is a complex character which is not afforded by a real representation, and
 3. `χ` is quaternionic character, i.e. it is real valued, but is not afforded a real representation.

In cases 2. and 3. `2re(χ) = χ + conj(χ)` corresponds to an irreducible character
afforded by a real representation.
"""
function frobenius_schur_indicator(χ::Character)
    @assert isirreducible(χ)

    pmap = powermap(χ)
    ι = sum(
        length(c) * χ[pmap[i, 2]] for (i, c) in enumerate(conjugacy_classes(χ))
    )

    ι_int = Int(ι)
    ordG = sum(length, conjugacy_classes(χ))
    d, r = divrem(ι_int, ordG)
    @assert r == 0 "Non integral Frobenius Schur Indicator: $(ι_int) = $d * $ordG + $r"
    return d
end

Base.isreal(χ::Character) = frobenius_schur_indicator(χ) > 0

function Base.show(io::IO, ::MIME"text/plain", χ::Character)
    println(io, "Character over ", eltype(χ))
    _print_char(io, χ)
end

Base.show(io::IO, χ::Character) = _print_char(io, χ)

function _print_char(io::IO, χ::Character)
    first = true
    for (i, c) in enumerate(constituents(χ))
        iszero(c) && continue
        first || print(io, " ")
        print(io, ((c < 0 || first) ? "" : '+'))
        !isone(c) && print(io, c, '·')
        print(io, 'χ', FiniteFields.subscriptify(i))
        first = false
    end
end

## characters defined by actions/homomorphisms

function _action_class_fun(
    conjugacy_cls::AbstractVector{CCl},
) where {CCl <: AbstractOrbit{<:AbstractPerm}}
    vals = Int[PermutationGroups.nfixedpoints(first(cc)) for cc in conjugacy_cls]
    # in general:
    # vals = [tr(matrix_representative(first(cc))) for cc in conjugacy_cls]
    return ClassFunction(vals, conjugacy_cls)
end

function _action_class_fun(
    conjugacy_cls::AbstractVector{CCl},
) where {CCl <: AbstractOrbit{<:AbstractMatrix}}
    vals = [tr(first(cc)) for cc in conjugacy_cls]
    return ClassFunction(vals, conjugacy_cls)
end

function _action_class_fun(
    hom::InducedActionHomomorphism{<:ByPermutations},
    conjugacy_cls
) where {CCl <: AbstractOrbit{<:AbstractPerm}}
    vals = Int[PermutationGroups.nfixedpoints(hom(first(cc))) for cc in conjugacy_cls]
    # in general:
    # vals = [tr(matrix_representative(first(cc))) for cc in conjugacy_cls]
    return ClassFunction(vals, conjugacy_cls)
end

function _action_class_fun(
    hom::InducedActionHomomorphism{<:ByLinearTransformation},
    conjugacy_cls,
)
    vals = [tr(hom(first(cc))) for cc in conjugacy_cls]
    return ClassFunction(vals, conjugacy_cls)
end

"""
    action_character(conjugacy_cls, tbl::CharacterTable)
Return the character of the representation given by the elements in the
conjugacy classes `conjugacy_cls`.

This corresponds to the classical definition of characters as a traces of the
representation matrices.
"""
function action_character(conjugacy_cls, tbl::CharacterTable)
    ac_char = _action_class_fun(conjugacy_cls)
    constituents = Int[dot(ac_char, χ) for χ in irreducible_characters(tbl)]
    return Character(tbl, constituents)
end

"""
    action_character(hom::InducedActionHomomorphism, tbl::CharacterTable)
Return the character of the representation given by the images under `hom` of
elements in `conjugacy_classes(tbl)`.

This corresponds to the classical definition of characters as a traces of the
representation matrices.
"""
function action_character(hom::InducedActionHomomorphism, tbl::CharacterTable)
    ac_char = _action_class_fun(hom, conjugacy_classes(tbl))
    constituents = Int[dot(ac_char, χ) for χ in irreducible_characters(tbl)]
    return Character(tbl, constituents)
end
