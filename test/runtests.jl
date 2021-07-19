using Test

using LinearAlgebra
using SymbolicWedderburn
using SymbolicWedderburn.FiniteFields
using LinearAlgebra

using GroupsCore
using PermutationGroups
using Cyclotomics

include("smallgroups.jl")

include("gf.jl")
include("eigenspacedecomposition.jl")
include("ccmatrix.jl")
include("dixon.jl")
include("characters.jl")
include("projections.jl")
include("action_permutation.jl")
include("action_linear.jl")

if VERSION >= v"1.6.0"
    @testset "Examples" begin
        using Pkg
        Pkg.activate(joinpath(@__DIR__, "..", "examples"))
        Pkg.instantiate()
        include("../examples/ex_C2_linear.jl")
        include("../examples/ex_C4.jl")
        include("../examples/ex_motzkin.jl")
    end
end
