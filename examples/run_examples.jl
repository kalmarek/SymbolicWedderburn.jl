using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()
include(joinpath(@__DIR__, "ex_C2_linear.jl"))
include(joinpath(@__DIR__, "ex_S4.jl"))
include(joinpath(@__DIR__, "ex_motzkin.jl"))
include(joinpath(@__DIR__, "ex_robinson_form.jl"))
