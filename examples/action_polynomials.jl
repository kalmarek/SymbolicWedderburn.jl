struct PermutingVariables <: SymbolicWedderburn.ByPermutations end

function SymbolicWedderburn.action(::PermutingVariables, g::PermutationGroups.AbstractPerm, m::Monomial)
    v = variables(m)
    return m(v => map(i -> v[i^g], eachindex(v)))
end

function SymbolicWedderburn.action(::PermutingVariables, g::PermutationGroups.AbstractPerm, v::AbstractVector)
    return map(i -> v[i^g], eachindex(v))
end
