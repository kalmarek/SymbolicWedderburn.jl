using Test

using LinearAlgebra
using SymbolicWedderburn
using SymbolicWedderburn.FiniteFields
using LinearAlgebra

using GroupsCore
using PermutationGroups
using Cyclotomics

include("smallgroups.jl")

@testset "Characters" begin
    import SymbolicWedderburn.Characters
    include("gf.jl")
    include("eigenspacedecomposition.jl")
    include("ccmatrix.jl")
    include("dixon.jl")
    include("characters.jl")
end
include("projections.jl")
include("sa_basis.jl")
include("action_permutation.jl")
include("action_linear.jl")
include("action_dihedral.jl")

if VERSION >= v"1.7.0" && !haskey(ENV, "CI")
    @testset "Examples" begin
        include("../examples/run_examples.jl")
    end
end
