abstract type AbstractAction{T} end

struct OnPoints{T} <: AbstractAction{T}
    features::Vector{T}
    reversef::Dict{T, Int}
end

Base.getindex(act::OnPoints, i::Integer) = act.features[i]
Base.getindex(act::OnPoints{T}, f::T) where T = act.reversef[f]
(act::OnPoints)(p::Perm{I}) where I = Perm(I[act[f^p] for f in act.features])

function (act::OnPoints)(orb::O) where {T, O<:AbstractOrbit{T, Nothing}}
    elts = act.(orb)
    dict = Dict(e=>nothing for e in elts)
    return O(elts, dict)
end
