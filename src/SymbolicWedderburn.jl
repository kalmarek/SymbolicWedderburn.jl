module SymbolicWedderburn

using LinearAlgebra
import AbstractAlgebra
using Primes

using PermutationGroups

include("gf.jl")
include("eigenspacedecomposition.jl")
include("ccmatrix.jl")
include("dixon.jl")

end # module
