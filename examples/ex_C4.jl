using SymbolicWedderburn
using PermutationGroups
using Cyclotomics

using SparseArrays

using DynamicPolynomials
using SumOfSquares
using SCS

OPTIMIZER = optimizer_with_attributes(
    SCS.Optimizer,
    "acceleration_lookback" => 10,
    "max_iters" => 10_000,
    "eps" => 1e-6,
    "linear_solver" => SCS.DirectSolver,
)

include(joinpath(@__DIR__, "action_polynomials.jl"))

N = 4

@polyvar x[1:N]

f = sum(x) +
    sum(x .^ 2) +
    (sum((x .+ 1) .^ 2 .* (x .+ x') .^ 2))^2 * (1 + sum(x .^ 2))

basis = monomials(x, 0:(DynamicPolynomials.maxdegree(f)÷2))

ts, obj, st = let f = f, basis = basis, m = SOSModel(OPTIMIZER)
    @variable m t
    @objective m Max t
    @variable m sos SOSPoly(basis)
    @constraint m f - t == sos
    optimize!(m)
    @info (m,) termination_status(m) objective_value(m) solve_time(m)
    termination_status(m), objective_value(m), solve_time(m)
end

ts_sa, obj_sa, st_sa = let f = f,
    basis = basis,
    m = SOSModel(OPTIMIZER),
    # G = PermGroup(Perm([2:N; 1])),
    G = PermGroup([perm"(1,2)", Perm([2:N; 1])])

    t = @timed let
        ssimple_basis = SymbolicWedderburn.symmetry_adapted_basis(Float64, G, basis, PermutingVariables())
        SparseMatrixCSC{Float64,Int}.(SymbolicWedderburn.basis.(ssimple_basis))
    end

    sa_basis, symmetry_adaptation_time = t.value, t.time

    let m = m, basis = basis, R = sa_basis
        @variable m t
        @objective m Max t

        soses = @variable m [r in R] SOSPoly(FixedPolynomialBasis(r * basis))
        @constraint m f - t == sum(soses)

        optimize!(m)
        @info (m,) termination_status(m) objective_value(m) solve_time(m) symmetry_adaptation_time
        termination_status(m), objective_value(m), solve_time(m)
    end
end

@assert ts == ts_sa == SumOfSquares.MOI.OPTIMAL
@assert isapprox(obj, obj_sa, atol=1e-4)
@assert st_sa < st
