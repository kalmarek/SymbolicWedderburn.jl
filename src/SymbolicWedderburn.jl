module SymbolicWedderburn

using LinearAlgebra
using Primes

using Cyclotomics
using GroupsCore
using StarAlgebras

import PermutationGroups:
    AbstractOrbit, AbstractPerm, AbstractPermutationGroup, Orbit, Perm, degree
import PermutationGroups

export conjugacy_classes, symmetry_adapted_basis

include("gf.jl")
include("eigenspacedecomposition.jl")
include("cmmatrix.jl")

include("powermap.jl")
include("character_tables.jl")
include("characters.jl")
include("dixon.jl")

include("actions.jl")
include("projections.jl")
include("minimal_projections.jl")
include("sa_basis.jl")

end # module
