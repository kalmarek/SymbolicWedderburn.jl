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
    memoize = false,
)
    hom = ExtensionHomomorphism(action, basis)
    cache = Dict(s => induce(hom, s) for s in gens(G))
    cache[one(G)] = induce(hom, one(G))
    sehom = SchreierExtensionHomomorphism(
        hom,
        cache,
        PG.SchreierTransversal(one(G), gens(G), *),
        memoize,
        Threads.SpinLock(),
    )

    return sehom
end
end
