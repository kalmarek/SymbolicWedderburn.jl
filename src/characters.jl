abstract type AbstractClassFunction{T} end # <: AbstractVector{T} ??

Base.eltype(::AbstractClassFunction{T}) where {T} = T

function LinearAlgebra.dot(
    χ::AbstractClassFunction{T},
    ψ::AbstractClassFunction{T},
) where {T}

    val = sum(
        length(cc) * χ[i] * ψ[-i]
        for (i, cc) in enumerate(conjugacy_classes(χ))
    )
    orderG = sum(length, conjugacy_classes(χ))
    val = div(val, orderG)
    return val
end

Base.isreal(χ::AbstractClassFunction) = all(isreal, values(χ))

####################################
# Characters

struct VirtualCharacter{T,CCl<:AbstractOrbit} <: AbstractClassFunction{T}
    vals::Vector{T}
    inv_of::Vector{Int}
    cc::Vector{CCl}
end

struct Character{T,CCl<:AbstractOrbit} <: AbstractClassFunction{T}
    vals::Vector{T}
    inv_of::Vector{Int}
    cc::Vector{CCl}
end

const ClassFunction = Union{Character,VirtualCharacter}

if VERSION >= v"1.3.0"
    function (χ::AbstractClassFunction)(g::PermutationGroups.GroupElem)
        for (i, cc) in enumerate(conjugacy_classes(χ))
            g ∈ cc && return χ[i]
        end
        throw(DomainError(
            g,
            "element does not belong to conjugacy classes of χ",
        ))
    end
else
    function (χ::Character)(g::PermutationGroups.GroupElem)
        for (i, cc) in enumerate(conjugacy_classes(χ))
            g ∈ cc && return χ[i]
        end
        throw(DomainError(
            g,
            "element does not belong to conjugacy classes of χ",
        ))
    end
    function (χ::VirtualCharacter)(g::PermutationGroups.GroupElem)
        for (i, cc) in enumerate(conjugacy_classes(χ))
            g ∈ cc && return χ[i]
        end
        throw(DomainError(
            g,
            "element does not belong to conjugacy classes of χ",
        ))
    end
end

function Character(
    vals::AbstractVector{T},
    ccls::AbstractVector{CCl},
) where {T,CCl<:AbstractOrbit}
    χ = Character{T,CCl}(vals, _inv_of(ccls), ccls)
    return χ
end

function VirtualCharacter(
    vals::AbstractVector{T},
    ccls::AbstractVector{CCl},
) where {T,CCl<:AbstractOrbit}
    χ = VirtualCharacter{T,CCl}(vals, _inv_of(ccls), ccls)
    return χ
end

VirtualCharacter(χ::Character{T,Cl}) where {T,Cl} = VirtualCharacter{T,Cl}(χ)
VirtualCharacter{T,Cl}(χ::Character{S,Cl}) where {T,S,Cl} =
    VirtualCharacter{T,Cl}(values(χ), χ.inv_of, conjugacy_classes(χ))

Base.values(χ::ClassFunction) = χ.vals
PermutationGroups.conjugacy_classes(χ::ClassFunction) = χ.cc

PermutationGroups.degree(χ::Character) =
    Int(χ(one(first(first(conjugacy_classes(χ))))))

PermutationGroups.degree(χ::VirtualCharacter) =
    χ(one(first(first(conjugacy_classes(χ)))))

Base.conj(χ::Cf) where {Cf<:ClassFunction} =
    Cf(values(χ)[χ.inv_of], χ.inv_of, conjugacy_classes(χ))

Base.@propagate_inbounds function Base.getindex(χ::ClassFunction, i::Integer)
    @boundscheck checkbounds(values(χ), abs(i))
    if i < 0
        return @inbounds values(χ)[χ.inv_of[abs(i)]]
    else
        return @inbounds values(χ)[i]
    end
end

Base.eltype(χ::ClassFunction) = eltype(values(χ))

function _inv_of(cc::AbstractVector{<:AbstractOrbit})
    inv_of = zeros(Int, size(cc))
    for (i, c) in enumerate(cc)
        g = inv(first(c))
        inv_of[i] = something(findfirst(k -> g in k, cc), 0)
    end
    any(iszero, inv_of) &&
        throw(ArgumentError("Could not find the conjugacy class for inverse of $(first(cc[findfirst(iszero, inv_of)]))."))
    return inv_of
end

function normalize!(χ::Character)

    ccG = conjugacy_classes(χ)
    id = one(first(first(ccG)))

    k = χ(id)
    if !isone(k)
        values(χ) .*= inv(k)
    end

    # ⟨χ, χ⟩ = 1/d²

    deg = sqrt(inv(dot(χ, χ)))
    # @debug "normalizing with" n dot(χ, χ) χ(id) χ

    # normalizing χ
    values(χ) .*= deg

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

            print(io, string($C) * " : ")
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
