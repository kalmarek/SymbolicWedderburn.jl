module SymbolicWedderburn

using LinearAlgebra
using SparseArrays
using Primes

using Cyclotomics
using GroupsCore
using PermutationGroups
using StarAlgebras

export symmetry_adapted_basis

macro spawn_compat(expr)
    @static if VERSION < v"1.3.0"
        return :(@async $(esc(expr)))
    else
        return :(Threads.@spawn $(esc(expr)))
    end
end

include("Characters/Characters.jl")
using .Characters
import .Characters: row_echelon_form!
import .Characters.FiniteFields

include("actions.jl")
include("action_characters.jl")
include("projections.jl")
include("minimal_projections.jl")
include("sa_basis.jl")

end # module
