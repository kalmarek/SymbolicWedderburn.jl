"""
    AbstractClassFunction
Abstract type representing functions constant on conjugacy classes of a group (e.g. characters).

Subtypes need to implement:
 * `getindex(χ, i::Integer)`: return the value on `i`-th conjugacy class.
 * Indexing with negative integers should return values on the class which contains inverses of the `i`-th class
 * `χ(g::GroupElem)`: return value of `χ` on a group element (only for julia < 1.3)
 * `values(χ)`: return an iterator over values
 * `conjugacy_classes(χ)`: return iterator over conjugacy classes of `χ`.
 * `conj(χ)`: return the conjugate function

It is assumed that two class functions on the same group will return **identical** (ie. `===`) conjugacy_classes.
"""
abstract type AbstractClassFunction{T} end # <: AbstractVector{T} ??

Base.eltype(::Type{<:AbstractClassFunction{T}}) where {T} = T

function LinearAlgebra.dot(χ::AbstractClassFunction, ψ::AbstractClassFunction)
    val = sum(
        length(cc) * χ[i] * ψ[-i] for
        (i, cc) in enumerate(conjugacy_classes(χ))
    )
    orderG = sum(length, conjugacy_classes(χ))
    val = div(val, orderG)
    return val
end

Base.:(==)(χ::AbstractClassFunction, ψ::AbstractClassFunction) =
    conjugacy_classes(χ) === conjugacy_classes(ψ) && values(χ) == values(ψ)

Base.hash(χ::AbstractClassFunction, h::UInt) =
    hash(conjugacy_classes(χ), hash(values(χ), hash(AbstractClassFunction, h)))

####################################
# Characters

struct VirtualCharacter{T,CCl} <: AbstractClassFunction{T}
    vals::Vector{T}
    inv_of::Vector{Int}
    cc::Vector{CCl}
end

struct Character{T,CCl} <: AbstractClassFunction{T}
    vals::Vector{T}
    inv_of::Vector{Int}
    cc::Vector{CCl}
end

const ClassFunction = Union{Character,VirtualCharacter}

"""
    affordable_real!(χ::ClassFunction[, pmap::PowerMap])
Return either `χ` or `2re(χ)` depending whether `χ` is afforded by a real representation, modifying `χ` in place.
"""
function affordable_real!(
    χ::ClassFunction,
    pmap = PowerMap(conjugacy_classes(χ)),
)
    ι = frobenius_schur_indicator(χ, pmap)
    if ι <= 0 # i.e. χ is complex or quaternionic
        for i in eachindex(χ.vals)
            χ.vals[i] += conj(χ.vals[i])
        end
    end
    return χ
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
function frobenius_schur_indicator(
    χ::AbstractClassFunction,
    pmap::PowerMap = PowerMap(conjugacy_classes(χ)),
)
    ι = sum(
        length(c) * χ[pmap[i, 2]] for (i, c) in enumerate(conjugacy_classes(χ))
    )

    ι_int = Int(ι)
    ordG = sum(length, conjugacy_classes(χ))
    d, r = divrem(ι_int, ordG)
    @assert r == 0 "Non integral Frobenius Schur Indicator: $(ι_int) = $d * $ordG + $r"
    return d
end

Base.isreal(χ::AbstractClassFunction) = frobenius_schur_indicator(χ) > 0

if VERSION >= v"1.3.0"
    function (χ::AbstractClassFunction)(g::GroupsCore.GroupElement)
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
end

###################################
# Specific definitions for ClassFunction
const ClassFunction = Union{Character,VirtualCharacter}

function Character(
    vals::AbstractVector{T},
    ccls::AbstractVector{CCl},
) where {T,CCl}
    χ = Character{T,CCl}(vals, _inv_of(ccls), ccls)
    return χ
end

function VirtualCharacter(
    vals::AbstractVector{T},
    ccls::AbstractVector{CCl},
) where {T,CCl}
    χ = VirtualCharacter{T,CCl}(vals, _inv_of(ccls), ccls)
    return χ
end

function Base.deepcopy_internal(χ::CF, dict::IdDict) where {CF<:ClassFunction}
    return CF(copy(χ.vals), χ.inv_of, conjugacy_classes(χ))
end

Character{S}(χ::Character{T,Cl}) where {S,T,Cl} =
    Character{S, Cl}(values(χ), χ.inv_of, conjugacy_classes(χ))

VirtualCharacter(χ::Character{T,Cl}) where {T,Cl} = VirtualCharacter{T,Cl}(χ)
VirtualCharacter{T}(χ::Character{S,Cl}) where {T,S,Cl} =
    VirtualCharacter{T,Cl}(χ)
VirtualCharacter(χ::VirtualCharacter) = deepcopy(χ)

VirtualCharacter{T,Cl}(χ::ClassFunction) where {T,S,Cl} =
    VirtualCharacter{T,Cl}(values(χ), χ.inv_of, conjugacy_classes(χ))

Base.@propagate_inbounds function Base.getindex(χ::ClassFunction, i::Integer)
    @boundscheck checkbounds(values(χ), abs(i))
    if i < 0
        return @inbounds values(χ)[χ.inv_of[abs(i)]]
    else
        return @inbounds values(χ)[i]
    end
end

if VERSION < v"1.3.0"
    function (χ::Character)(g::GroupsCore.GroupElement)
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
    function (χ::VirtualCharacter)(g::GroupsCore.GroupElement)
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
end

Base.values(χ::ClassFunction) = χ.vals
conjugacy_classes(χ::ClassFunction) = χ.cc

Base.conj(χ::Cf) where {Cf<:ClassFunction} =
    Cf(values(χ)[χ.inv_of], χ.inv_of, conjugacy_classes(χ))

PermutationGroups.degree(χ::Character) = Int(χ.vals[1])
    # Int(χ(one(first(first(conjugacy_classes(χ))))))

function _inv_of(cc::AbstractVector{<:AbstractOrbit})
    inv_of = zeros(Int, size(cc))
    for (i, c) in enumerate(cc)
        g = inv(first(c))
        inv_of[i] = something(findfirst(k -> g in k, cc), 0)
    end
    any(iszero, inv_of) && throw(
        ArgumentError(
            "Could not find the conjugacy class for inverse of $(first(cc[findfirst(iszero, inv_of)])).",
        ),
    )
    return inv_of
end

function normalize!(χ::Character)

    ccG = conjugacy_classes(χ)
    id = one(first(first(ccG)))

    k = χ(id)
    if !isone(k)
        χ.vals .*= inv(k)
    end

    # ⟨χ, χ⟩ = 1/d²

    deg = sqrt(inv(dot(χ, χ)))
    @debug "normalizing with" n dot(χ, χ) χ(id) χ

    # normalizing χ
    χ.vals .*= deg

    return χ
end

for C in (:Character, :VirtualCharacter)
    @eval begin
        function Base.show(io::IO, ::MIME"text/plain", χ::$C)
            println(io, string($C) * " over $(eltype(χ))")
            cc_reps = string.(first.(conjugacy_classes(χ)))
            k = maximum(length.(cc_reps))

            vals = values(χ)
            for i = 1:length(vals)
                (c, v) = cc_reps[i], vals[i]
                if i == length(values(χ))
                    print(io, rpad("$c^G", k + 3), "→ \t", v)
                else
                    print(io, rpad("$c^G", k + 3), "→ \t", v, "\n")
                end
            end
        end

        function Base.show(io::IO, χ::$C)
            v = values(χ)
            io = IOContext(io, :typeinfo => eltype(v))
            limited = get(io, :limit, false)
            opn, cls = '[', ']'

            print(io, string($C) * ": ")
            if limited && length(v) > 20
                Base.show_delim_array(io, v, opn, ",", "", false, 1, 9)
                print(io, "  …  ")
                l = length(v)
                Base.show_delim_array(io, v, "", ",", cls, false, l - 9, l)
            else
                Base.show_delim_array(io, v, opn, ",", cls, false)
            end
        end
    end
end
