module SymbolicWedderburn

using LinearAlgebra
import AbstractAlgebra
using Primes

using Pkg
Pkg.add(PackageSpec(url="https://github.com/kalmarek/PermutationGroups.jl"))
using PermutationGroups

include("gf.jl")
include("eigenspacedecomposition.jl")
include("ccmatrix.jl")

include("characters.jl")
include("powermap.jl")
include("dixon.jl")

end # module
