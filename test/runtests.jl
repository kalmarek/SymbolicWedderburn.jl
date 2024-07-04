using Test

using LinearAlgebra
using SparseArrays

using GroupsCore
import AbstractPermutations as AP
import PermutationGroups as PG
using Cyclotomics

using SymbolicWedderburn
import SymbolicWedderburn as SW
using SymbolicWedderburn.FiniteFields
import SymbolicWedderburn.SA as SA
using SymbolicWedderburn.StarAlgebras

# action on free words
include("free_words.jl")
include("action_permutation.jl")

# actions on polynomials
include(joinpath(dirname(@__DIR__), "examples", "action_polynomials.jl"))
include(joinpath(dirname(@__DIR__), "examples", "dihedral.jl"))
include(joinpath(dirname(@__DIR__), "examples", "sos_problem.jl"))

include("action_linear.jl")
include("action_dihedral.jl")
include("action_invalid.jl")

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

if VERSION >= v"1.7.0" && !haskey(ENV, "CI")
    @testset "Examples" begin
        include("../examples/run_examples.jl")
    end
end

