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
include(joinpath(dirname(@__DIR__), "examples", "solver.jl"))

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

# Small groups generation:
#= GAP code:
H := [];
for i in [1..63] do
    Add(H, List(AllSmallGroups(i), G->Image(IsomorphismPermGroup(G))));
od;
PrintTo("/tmp/groups.gap", H);
=#

#=julia code
GAPgroups_str = join(readlines("/tmp/groups.gap"), "");
GAPgroups_str = replace(GAPgroups_str, "Group"=>"\nPermGroup");
GAPgroups_str = replace(GAPgroups_str, r" *"=>"");
perm_regex = r"((\(\d+(,\d+)*\)?)+)";
let fn = joinpath(@__DIR__, "smallgroups.jl")
    open(fn, "w") do file

        print(file, """
        import PermutationGroups: PermGroup, @perm_str

        const SmallPermGroups = """)
        println(file, replace(GAPgroups_str, perm_regex=> s"perm\"\1\""))
    end
    read(fn, String)
end
=#

