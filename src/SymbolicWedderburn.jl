module SymbolicWedderburn

using LinearAlgebra
using Primes

using Cyclotomics
using GroupsCore

import PermutationGroups:
    AbstractOrbit, AbstractPerm, AbstractPermutationGroup, Orbit, Perm, degree
import PermutationGroups

export conjugacy_classes, symmetry_adapted_basis

include("gf.jl")
include("eigenspacedecomposition.jl")
include("cmmatrix.jl")

include("powermap.jl")
include("characters.jl")
include("characters_arith.jl")
include("dixon.jl")

include("actions.jl")
include("projections.jl")

end # module
