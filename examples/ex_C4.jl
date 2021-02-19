using SymbolicWedderburn
using PermutationGroups
using Cyclotomics

using SparseArrays

using DynamicPolynomials
using SumOfSquares
using SCS

include(joinpath(@__DIR__, "action_polynomials.jl"))

SCS_Indirect, SCS_Direct =
    let params = (
            "acceleration_lookback" => 10,
            "max_iters" => 10_000,
            "eps" => 1e-6,
        )
        indir = optimizer_with_attributes(
            SCS.Optimizer,
            params...,
            "linear_solver" => SCS.IndirectSolver,
        )

        dir = optimizer_with_attributes(
            SCS.Optimizer,
            params...,
            "linear_solver" => SCS.DirectSolver,
        )
        indir, dir
    end

OPTIMIZER = let
    SCS_Direct
end

N = 4

@polyvar x[1:N]

f =
    sum(x) +
    sum(x .^ 2) +
    (sum((x .+ 1) .^ 2 .* (x .+ x') .^ 2))^2 * (1 + sum(x .^ 2))

basis = monomials(x, 0:(DynamicPolynomials.maxdegree(f)÷2))

@time let f = f, basis = basis, m = SOSModel(OPTIMIZER)
    @variable m t
    @objective m Max t
    @variable m sos SOSPoly(basis)
    @constraint m f - t == sos
    optimize!(m)
    @info (m,) termination_status(m) objective_value(m) solve_time(m)
end

@time let f = f,
    basis = basis,
    m = SOSModel(OPTIMIZER),
    # G = PermGroup(Perm([2:N; 1])),
    G = PermGroup([perm"(1,2)", Perm([2:N; 1])])

    t = @timed let
        sa_basis = SymbolicWedderburn.symmetry_adapted_basis(G, exponents.(basis), permuted)
        SparseMatrixCSC{Float64,Int}.(sa_basis)
    end

    sa_basis, symmetry_adaptation_time = t.value, t.time

    let m = m, basis = basis, R = sa_basis
        @variable m t
        @objective m Max t

        soses = @variable m [r in R] SOSPoly(FixedPolynomialBasis(r * basis))
        @constraint m f - t == sum(soses)

        optimize!(m)
        @info (m,) termination_status(m) objective_value(m) solve_time(m) symmetry_adaptation_time
    end
end
