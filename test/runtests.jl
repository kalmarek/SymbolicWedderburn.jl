using Test

using SymbolicWedderburn
using SymbolicWedderburn.FiniteFields

using Pkg
Pkg.add(PackageSpec(url="https://github.com/kalmarek/PermutationGroups.jl"))
using PermutationGroups


include("gf.jl")
include("eigenspacedecomposition.jl")
include("ccmatrix.jl")
include("dixon.jl")
