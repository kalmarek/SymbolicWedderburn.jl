using SymbolicWedderburn
using PermutationGroups
using Cyclotomics

using SparseArrays

using DynamicPolynomials
using SumOfSquares
using SCS

const OPTIMIZER = optimizer_with_attributes(
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

no_symmetry = let f = f, basis = monomials(x, 0:(DynamicPolynomials.maxdegree(f)÷2))
    stats = @timed let m = SOSModel(OPTIMIZER)
        @variable m t
        @objective m Max t
        @variable m sos SOSPoly(basis)
        @constraint m f - t == sos
        m
    end
    m = stats.value
    optimize!(m)
    @info (m,) termination_status(m) objective_value(m) creation_time=(stats.time) solve_time(m)
    (
        status=termination_status(m),
        objective=objective_value(m),
        creation_t=stats.time,
        solve_t=solve_time(m)
    )
end

semisimple_dec = let f = f, G = PermGroup([perm"(1,2)", Perm([2:N; 1])])
    basis = monomials(x, 0:(DynamicPolynomials.maxdegree(f)÷2))

    stats = @timed begin
        sa_basis, symmetry_adaptation_time, = @timed let
            ssimple_basis = SymbolicWedderburn.symmetry_adapted_basis(
                Rational{Int},
                G,
                VariablePermutation(),
                basis,
                semisimple=true
            )
            sp = sparse.(ssimple_basis)
            droptol!.(sp, 1e-12)
        end

        m = SOSModel(OPTIMIZER)

        let basis = basis, R = sa_basis
            @variable m t
            @objective m Max t

            soses = @variable m [r in R] SOSPoly(FixedPolynomialBasis(r * basis))
            @constraint m f - t == sum(soses)
        end
        m
    end
    m = stats.value
    optimize!(m)
    @info (m,) termination_status(m) objective_value(m) solve_time(m) creation_time=(stats.time) symmetry_adaptation_time
    (
        status=termination_status(m),
        objective=objective_value(m),
        creation_t=stats.time,
        solve_t=solve_time(m)
    )
end

@assert no_symmetry.status == MOI.OPTIMAL
@assert semisimple_dec.status == MOI.OPTIMAL
@assert isapprox(no_symmetry.objective, semisimple_dec.objective, atol=1e-4)
@assert semisimple_dec.solve_t < no_symmetry.solve_t

# This is faster, but not quite the speedup we'd like to see.
# Let us try now to use the decomposition into simple summands
include(joinpath(@__DIR__, "util.jl"))

wedderburn_dec = let f=f, G = PermGroup([perm"(1,2)", Perm([2:N; 1])]), T = Float64

    basis_half = monomials(x, 0:(DynamicPolynomials.maxdegree(f)÷2))

    wedderburn, symmetry_adaptation_time, = @timed WedderburnDecomposition(
        T,
        G,
        VariablePermutation(),
        monomials(x, 0:DynamicPolynomials.maxdegree(f)),
        basis_half
    )

    M = [SymbolicWedderburn.basis(wedderburn)[x*y] for x in basis_half, y in basis_half]

    stats = @timed begin
        m = let m = JuMP.Model(OPTIMIZER)
            JuMP.@variable m t
            JuMP.@objective m Max t
            psds = [
                JuMP.@variable(m, [1:d, 1:d] in PSDCone())
                for d in size.(direct_summands(wedderburn), 1)
            ]

            # preallocating
            Mπs = zeros.(T, size.(psds))
            M_orb = similar(M, T)
            tmps = SymbolicWedderburn._tmps(wedderburn)

            C = DynamicPolynomials.coefficients(f-t, SymbolicWedderburn.basis(wedderburn))

            for iv in invariant_vectors(wedderburn)
                c = dot(C, iv)
                # C needs to be invariant under G,
                # which is eqivalent to being constant on the orbits
                let cc = C[findfirst(!iszero, iv)]
                    @assert all(C[m] == cc for m in findall(!iszero, iv))
                end

                # invariant_constraint! is defined locally in util.jl
                # it essentially amounts to averaging constraint matrices defined by M
                # with weights provided by iv and storing the result in (dense) matrix M_orb
                M_orb = invariant_constraint!(M_orb, M, iv)
                Mπs = SymbolicWedderburn.diagonalize!(Mπs, M_orb, wedderburn, tmps)

                JuMP.@constraint m sum(
                    dot(Mπ, Pπ) for (Mπ, Pπ) in zip(Mπs, psds) if !iszero(Mπ)
                ) == c
            end
            m
        end
    end

    m = stats.value
    optimize!(m)
    @info (m,) termination_status(m) objective_value(m) solve_time(m) creation_time=(stats.time) symmetry_adaptation_time
    (
        status=termination_status(m),
        objective=objective_value(m),
        creation_t=stats.time,
        solve_t=solve_time(m)
    )
end

@assert wedderburn_dec.status == MOI.OPTIMAL
@assert isapprox(wedderburn_dec.objective, semisimple_dec.objective, atol=1e-4)
@assert semisimple_dec.solve_t > wedderburn_dec.solve_t
@assert wedderburn_dec.solve_t < 2.0
