using SymbolicWedderburn
using PermutationGroups

using SparseArrays

using DynamicPolynomials
using SumOfSquares
using SCS

include(joinpath(@__DIR__, "action_polynomials.jl"))

OPTIMIZER = optimizer_with_attributes(
    SCS.Optimizer,
    "acceleration_lookback" => 0,
    "max_iters" => 20_000,
    "eps" => 1e-10,
    "linear_solver" => SCS.DirectSolver,
)

@polyvar x y z

motzkin = x^4 * y^2 + y^4 * x^2 - 3 * x^2 * y^2 + 1
g = (x^2 + y^2 + 1)
basis = monomials([x, y], 0:7)

@time let f = motzkin, basis = basis, m = SOSModel(OPTIMIZER)
    @variable m t >=0
    # @objective m Max t
    @variable m sos SOSPoly(basis)
    @constraint m g*f - t == sos
    optimize!(m)
    @info (m,) termination_status(m) objective_value(m) solve_time(m)
end

@time let f = motzkin,
    basis = basis,
    m = SOSModel(OPTIMIZER),
    G = PermGroup(perm"(1,2)")

    t = @timed let
        ssimple_basis = SymbolicWedderburn.symmetry_adapted_basis(G, exponents.(basis), PermutingVariables())
        SparseMatrixCSC{Float64,Int}.(SymbolicWedderburn.basis.(ssimple_basis))
    end

    sa_basis, symmetry_adaptation_time = t.value, t.time

    let m = m, basis = basis, R = sa_basis
        @variable m t >= 0
        # @objective m Max t

        soses = @variable m [r in R] SOSPoly(FixedPolynomialBasis(r * basis))
        @constraint m f - t == sum(soses)

        optimize!(m)
        @info m, termination_status(m) objective_value(m) solve_time(m) symmetry_adaptation_time
    end
end
