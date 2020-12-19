module SymbolicWedderburn

using LinearAlgebra
using Primes

using PermutationGroups
using Cyclotomics

include("gf.jl")
include("eigenspacedecomposition.jl")
include("cmmatrix.jl")

include("characters.jl")
include("characters_arith.jl")
include("powermap.jl")
include("dixon.jl")

include("actions.jl")
include("projections.jl")

end # module
