function __schreier_tree(G::Group)
    # stores (g => (h,s)) when g = h·s for a generator s of G;
    # this will be replaced by
    # `PermutationGroups.SchreierTransversal(one(G), gens(G), *)`
    # once it is fixed upstream

    id = one(G)
    orbit = [id]
    tree = Dict(id => (id, id))
    S = gens(G)
    for h in orbit
        for s in S
            g = h * s # I think we want action on the RIGHT here
            # so that reconstruction happens without reversal
            haskey(tree, g) && continue
            push!(orbit, g)
            tree[g] = (h, s)
        end
    end
    return tree
end

struct SchreierExtensionHomomorphism{
    A,
    T,
    E<:InducedActionHomomorphism{A,T},
    Ch,
    St,
} <: InducedActionHomomorphism{A,T}
    ehom::E
    cache::Ch
    schreier_tree::St
    memoize::Bool
    lock::Base.Threads.SpinLock
end

StarAlgebras.basis(h::SchreierExtensionHomomorphism) = basis(h.ehom)
action(h::SchreierExtensionHomomorphism) = action(h.ehom)

function SchreierExtensionHomomorphism(
    G::Group,
    action::Action,
    basis;
    memoize::Bool = false,
)
    return SchreierExtensionHomomorphism(
        _int_type(action),
        G,
        action,
        basis;
        memoize = memoize,
    )
end

function SchreierExtensionHomomorphism(
    ::Type{I},
    G::Group,
    action::Action,
    basis;
    memoize = false,
) where {I}
    hom = ExtensionHomomorphism(I, action, basis)
    cache = Dict(s => induce(hom, s) for s in gens(G))
    cache[one(G)] = induce(hom, one(G)) # needed?
    schr_hom = SchreierExtensionHomomorphism(
        hom,
        cache,
        __schreier_tree(G),
        memoize,
        Threads.SpinLock(),
    )

    return schr_hom
end
