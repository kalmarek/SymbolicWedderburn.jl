using LinearAlgebra
using SparseArrays

using GroupsCore
import SymbolicWedderburn as SW
import SymbolicWedderburn.StarAlgebras as SA
using PermutationGroups
import PermutationGroups.AP as AP

include(joinpath(@__DIR__, "action_polynomials.jl"))
include(joinpath(@__DIR__, "sos_problem.jl"))
include(joinpath(@__DIR__, "solver.jl"))
