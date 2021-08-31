module SymbolicWedderburn

using LinearAlgebra
using SparseArrays
using Primes

using Cyclotomics
using GroupsCore
using StarAlgebras

import PermutationGroups:
    AbstractOrbit, AbstractPerm, AbstractPermutationGroup, Orbit, Perm, degree
import PermutationGroups

export conjugacy_classes, symmetry_adapted_basis

macro spawn_compat(expr)
    @static if VERSION < v"1.3.0"
        return :(@async $(esc(expr)))
    else
        return :(Threads.@spawn $(esc(expr)))
    end
end

include("gf.jl")
include("eigenspacedecomposition.jl")
include("cmmatrix.jl")

include("actions.jl")
include("powermap.jl")
include("character_tables.jl")
include("characters.jl")
include("dixon.jl")

include("projections.jl")
include("minimal_projections.jl")
include("sa_basis.jl")

end # module
