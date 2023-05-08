abstract type Action end

abstract type InducedActionHomomorphism{A,T} end

#=
implements:
* basis(action_hom)
* action(hom::InducedActionHomomorphism) → Action
=#

Base.getindex(hom::InducedActionHomomorphism, i::Integer) = basis(hom)[i]
function Base.getindex(hom::InducedActionHomomorphism{A,T}, f::T) where {A,T}
    return basis(hom)[f]
end

coeff_type(hom::InducedActionHomomorphism) = coeff_type(action(hom))
_int_type(::Type{<:StarAlgebras.AbstractBasis{T,I}}) where {T,I} = I
_int_type(hom::InducedActionHomomorphism) = _int_type(typeof(basis(hom)))

# Exceeding typemax(UInt16) here would mean e.g. that you're trying to block-diagonalize
# an SDP constraint of size 65535×65535, which is highly unlikely ;)
_int_type(::Type{<:InducedActionHomomorphism}) = UInt16

function induce(hom::InducedActionHomomorphism, g::GroupElement)
    return induce(action(hom), hom, g)
end

function induce(ac::Action, hom::InducedActionHomomorphism, g::GroupElement)
    return throw(
        """No fallback is provided for $(typeof(ac)).
        You need to implement
        `induce(::$(typeof(ac)), ::$(typeof(hom)), ::$(typeof(g)))`.""",
    )
end

struct ExtensionHomomorphism{A<:Action,T,B<:StarAlgebras.AbstractBasis{T}} <:
       InducedActionHomomorphism{A,T}
    action::A
    basis::B
end

function ExtensionHomomorphism(action::Action, basis)
    return ExtensionHomomorphism(
        _int_type(ExtensionHomomorphism),
        action,
        basis,
    )
end

function ExtensionHomomorphism(
    ::Type{I},
    action::Action,
    basis,
) where {I<:Integer}
    return ExtensionHomomorphism(action, StarAlgebras.Basis{I}(basis))
end

# interface:
StarAlgebras.basis(hom::ExtensionHomomorphism) = hom.basis
action(hom::ExtensionHomomorphism) = hom.action

struct CachedExtensionHomomorphism{A,T,G,H,E<:InducedActionHomomorphism{A,T}} <:
       InducedActionHomomorphism{A,T}
    ehom::E
    cache::Dict{G,H}
    lock::Base.Threads.SpinLock
end

function CachedExtensionHomomorphism{G,H}(
    hom::InducedActionHomomorphism,
) where {G,H}
    return CachedExtensionHomomorphism(hom, Dict{G,H}(), Threads.SpinLock())
end

StarAlgebras.basis(h::CachedExtensionHomomorphism) = basis(h.ehom)
action(h::CachedExtensionHomomorphism) = action(h.ehom)

function CachedExtensionHomomorphism(
    G::Group,
    action::Action,
    basis;
    precompute = false,
)
    hom = ExtensionHomomorphism(action, basis)
    S = typeof(induce(hom, one(G)))
    chom = CachedExtensionHomomorphism{eltype(G),S}(hom)
    @sync if precompute
        for g in G
            Threads.@spawn begin
                induce(action, chom, g)
            end
        end
    end
    return chom
end

function induce(ac::Action, chom::CachedExtensionHomomorphism, g::GroupElement)
    return _induce(ac, chom, g)
end

function _induce(
    action::Action,
    chom::CachedExtensionHomomorphism,
    g::GroupElement,
)
    if !haskey(chom.cache, g)
        val = induce(action, chom.ehom, g)
        lock(chom.lock) do
            return chom.cache[g] = val
        end
    end
    return chom.cache[g]
end
