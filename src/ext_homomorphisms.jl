"""
    Top of the type hierarchy representing a group G acting on some object Ω.
    That is A:Ω×G->Ω or A:G×Ω->Ω (right/left actions) satisfying usual axioms.
    This induces the group homomorphism φ:G->Aut(Ω) called action homomorphism.
    SymbolicWedderburn is agnostic to right/left actions even though the method
    action(a::Action, g::GroupElement, ω) appears to imply a left action g⋅ω.
    Warning: no check on action(...) is performed!
""" 
abstract type Action end
abstract type InducedActionHomomorphism{A,T} end

Base.getindex(hom::InducedActionHomomorphism, i::Integer) = basis(hom)[i]
Base.getindex(hom::InducedActionHomomorphism{A,T}, f::T) where {A,T} =
    basis(hom)[f]

coeff_type(hom::InducedActionHomomorphism) = coeff_type(action(hom))
_int_type(::Type{<:StarAlgebras.AbstractBasis{T,I}}) where {T,I} = I
_int_type(hom::InducedActionHomomorphism) = _int_type(typeof(basis(hom)))

# Exceeding typemax(UInt16) here would mean e.g. that you're trying to block-diagonalize
# an SDP constraint of size 65535×65535, which is highly unlikely ;)
_int_type(::Type{<:InducedActionHomomorphism}) = UInt16

"""
    induce(hom::InducedActionHomomorphism, g::GroupElement)
Returns the induced action homomorphism at the group element g, i.e. φ(g).
"""
induce(hom::InducedActionHomomorphism, g::GroupElement) =
    induce(action(hom), hom, g)

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

ExtensionHomomorphism(action::Action, basis) =
    ExtensionHomomorphism(_int_type(ExtensionHomomorphism), action, basis)

ExtensionHomomorphism(::Type{I}, action::Action, basis) where {I<:Integer} =
    ExtensionHomomorphism(action, StarAlgebras.Basis{I}(basis))

# interface:
StarAlgebras.basis(hom::ExtensionHomomorphism) = hom.basis
action(hom::ExtensionHomomorphism) = hom.action

# Cache version: store the map φ in a dictionary with keys g ∈ G
struct CachedExtensionHomomorphism{A,T,G,H,E<:InducedActionHomomorphism{A,T}} <:
       InducedActionHomomorphism{A,T}
    ehom::E
    cache::Dict{G,H}
    lock::Base.Threads.SpinLock
end

CachedExtensionHomomorphism{G,H}(hom::InducedActionHomomorphism) where {G,H} =
    CachedExtensionHomomorphism(hom, Dict{G,H}(), Threads.SpinLock())

# interface:
StarAlgebras.basis(h::CachedExtensionHomomorphism) = basis(h.ehom)
action(h::CachedExtensionHomomorphism) = action(h.ehom)

# main constructor, default to lazy
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

induce(ac::Action, chom::CachedExtensionHomomorphism, g::GroupElement) =
    _induce(ac, chom, g)

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