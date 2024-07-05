using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

module Examples
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

include(joinpath(@__DIR__, "ex_C2_linear.jl"))
include(joinpath(@__DIR__, "ex_S4.jl"))
include(joinpath(@__DIR__, "ex_motzkin.jl"))
include(joinpath(@__DIR__, "ex_robinson_form.jl"))

end
