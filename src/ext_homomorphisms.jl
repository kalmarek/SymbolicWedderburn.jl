"""
    An abstract type representing a group G acting on some object Ω, aiding dispatch.
    that is a map A:Ω×G->Ω satisfying action axioms...somehow not quite right...
    The domain implies we work with RIGHT action, but instance of this type depends on user def
"""# should follow GAP and urge developer to use right action?
"""
    Top of the group action type hierarchy containing both left and right actions
"""
abstract type Action end # not parametric

"""
    An abstract type representing the induced homomorphism of a given Action
    that is the (group) homomorphism G->End(O)
"""
abstract type InducedActionHomomorphism{A,T} end # declares a family of parametric types
#=
    effectively? abstract type InducedActionHomomorphism end with all children
    abstract type InducedActionHomomorphism{A,T} end for all A and all T
    Do we declare InducedActionHomomorphism type?
=# 

"""
    Terminology from Hulpke: right action, notation exponential ω^g
    If G acts on set Ω, the action homomorphism is φ:G->SymmetricGroup(Ω), a group homomorphism
    In theory, two ways to specify: either give φ or give actfun below.
    In GAP, ω^g = actfun(ω, g) defines an action, i.e. actfun:Ω×G->G.
    In GAP, ActionHomomorphism returns a homomorphism G -> End(Ω), that is φ
"""


Base.getindex(hom::InducedActionHomomorphism, i::Integer) = basis(hom)[i] # why?
Base.getindex(hom::InducedActionHomomorphism{A,T}, f::T) where {A,T} =
    basis(hom)[f] # this covers L13

coeff_type(hom::InducedActionHomomorphism) = coeff_type(action(hom))
_int_type(::Type{<:StarAlgebras.AbstractBasis{T,I}}) where {T,I} = I
_int_type(hom::InducedActionHomomorphism) = _int_type(typeof(basis(hom)))

# Exceeding typemax(UInt16) here would mean e.g. that you're trying to block-diagonalize
# an SDP constraint of size 65535×65535, which is highly unlikely ;)
_int_type(::Type{<:InducedActionHomomorphism}) = UInt16

"""
    Computes the induced action homomorphism at the group element g, providing
    the image of G->End(Ω) at the point g, i.e. an endomorphism of Ω
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


"""
    Concrete type for actions on modules (implemented in Julia as a StarAlgebra).
    Here for each g ∈ G, this reduces to a map of basis of the StarAlgebra.
        ? extend action to basis ?
"""
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

"""
    Cache versions use a dictionary with keys g ∈ G and store all images of
    of the map G -> End(Ω) in the field cache
"""
struct CachedExtensionHomomorphism{A,T,G,H,E<:InducedActionHomomorphism{A,T}} <:
       InducedActionHomomorphism{A,T}
    ehom::E
    cache::Dict{G,H}
    lock::Base.Threads.SpinLock
end

CachedExtensionHomomorphism{G,H}(hom::InducedActionHomomorphism) where {G,H} =
    CachedExtensionHomomorphism(hom, Dict{G,H}(), Threads.SpinLock())

StarAlgebras.basis(h::CachedExtensionHomomorphism) = basis(h.ehom)
action(h::CachedExtensionHomomorphism) = action(h.ehom)

"""
    External constructor compute and cache G->End(M) of module M in ehom field
"""
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
