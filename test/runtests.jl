using Test

using LinearAlgebra
using SparseArrays

using GroupsCore
using PermutationGroups
using Cyclotomics

using SymbolicWedderburn
using SymbolicWedderburn.FiniteFields
using SymbolicWedderburn.StarAlgebras

include("action_permutation.jl")
include("action_linear.jl")
include("action_dihedral.jl")
include("action_invalid.jl")

if VERSION >= v"1.7.0" && !haskey(ENV, "CI")
    @testset "Examples" begin
        include("../examples/run_examples.jl")
    end
end

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
