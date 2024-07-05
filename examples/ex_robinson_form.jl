# # Dihedral symmetry of the Robinson form

# Based on version from SumOfSquares by Benoît Legat
# **Adapted from**: Example 5.4 of [GP04]
#
# [GP04] Gatermann, Karin and Parrilo, Pablo A.
# *Symmetry groups, semidefinite programs, and sums of squares*.
# Journal of Pure and Applied Algebra 192.1-3 (2004): 95-128.

using LinearAlgebra
using SparseArrays

using SymbolicWedderburn
import SymbolicWedderburn.SA as SA
using DynamicPolynomials

# definitions of general actions on polynomials:
include(joinpath(@__DIR__, "action_polynomials.jl"))
include(joinpath(@__DIR__, "sos_problem.jl"))
include(joinpath(@__DIR__, "solver.jl"))

@polyvar x y
const robinson_form =
    x^6 + y^6 - x^4 * y^2 - y^4 * x^2 - x^4 - y^4 - x^2 - y^2 + 3x^2 * y^2 + 1

# The Robinson form is invariant under the following action of the Dihedral group on monomials:
# The action of each element of the groups is to map the variables `x, y` to:
#
# | id | rotation | reflection |
# |----|----------|------------|
# | 0  | x, y     | y, x       |
# | 1  | -y, x    | -x, y      |
# | 2  | -x, -y   | -y, -x     |
# | 3  | y, -x    | x, -y      |

# implementation of the dihedral group:
include(joinpath(@__DIR__, "dihedral.jl"))

struct DihedralAction <: OnMonomials end

SymbolicWedderburn.coeff_type(::DihedralAction) = Int
function SymbolicWedderburn.action(
    ::DihedralAction,
    el::DihedralElement,
    mono::AbstractMonomial,
)
    x, y = variables(mono)
    if iseven(el.reflection + el.id)
        var_x, var_y = x, y
    else
        var_x, var_y = y, x
    end
    sign_x = 1 <= el.id <= 2 ? -1 : 1
    sign_y = 2 <= el.id ? -1 : 1
    return mono([x, y] => [sign_x * var_x, sign_y * var_y])
end

G = DihedralGroup(4)
for g in G
    @assert SymbolicWedderburn.action(DihedralAction(), g, robinson_form) ==
            robinson_form
end

OPTIMIZER =
    scs_optimizer(; max_iters = 10_000, alpha = 1.9, accel = -10, eps = 1e-5)

# A simple Gramm-matrix like formulation of SOS problem:
# * a single large pdd constraint and
# * one linear constraint for each monomial in the monomial basis of f.
no_symmetry = let f = robinson_form
    m, creation_t = @timed sos_problem(f)
    JuMP.set_optimizer(m, OPTIMIZER)
    optimize!(m)

    @info isapprox(objective_value(m), -3825 / 4096, rtol = 1e-4)

    (
        status = termination_status(m),
        objective = objective_value(m),
        symmetry_adaptation_t = 0.0,
        creation_t = creation_t,
        solve_t = solve_time(m),
    )
end

# Let's decompose basis into semi-simple summands and prepare the approximately
# reduced optimization problem. Since this is an isomorphism we reduce the
# number of variables from N² to n₁² + … + nₖ², where Σᵢ nᵢ = N. That is instead
# of one psd constraint of size N we have k constraints of sizes n₁,… nₖ.
# The number of linear constraints remains the same as before.

# Now given that f must be constant on orbits of monomials under the action of G
# let's exploit that to reduce the number of constraints (one for each orbit).
# Here we keep the single large psd constraint of size N×N.
orbit_dec = let f = robinson_form, T = Float64
    m, stats = sos_problem(
        f,
        DihedralGroup(4),
        DihedralAction();
        decompose_psd = false,
    )

    JuMP.set_optimizer(m, OPTIMIZER)
    optimize!(m)

    @info isapprox(objective_value(m), -3825 / 4096, rtol = 1e-4)

    (
        status = termination_status(m),
        objective = objective_value(m),
        symmetry_adaptation_t = stats["symmetry_adaptation"],
        creation_t = stats["model_creation"],
        solve_t = solve_time(m),
    )
end

semisimple_dec = let f = robinson_form, T = Float64
    m, stats = sos_problem(f, DihedralGroup(4), DihedralAction(); semisimple = true)
    JuMP.set_optimizer(m, OPTIMIZER)
    optimize!(m)

    @info "objective value/true objective: $(objective_value(m)) / $(-3825 / 4096)" isapprox(
        objective_value(m),
        -3825 / 4096,
        rtol = 1e-4,
    )

    (
        status = termination_status(m),
        objective = objective_value(m),
        symmetry_adaptation_t = stats["symmetry_adaptation"],
        creation_t = stats["model_creation"],
        solve_t = solve_time(m),
    )
end

wedderburn_dec = let f = robinson_form, T = Float64
    m, stats = sos_problem(
        f,
        DihedralGroup(4),
        DihedralAction();
        semisimple = false,
    )
    JuMP.set_optimizer(m, OPTIMIZER)
    optimize!(m)

    (
        status = termination_status(m),
        objective = objective_value(m),
        symmetry_adaptation_t = stats["symmetry_adaptation"],
        creation_t = stats["model_creation"],
        solve_t = solve_time(m),
    )
end

@info "Summary of Example: Robinson form" no_symmetry orbit_dec semisimple_dec wedderburn_dec

@assert wedderburn_dec.status == no_symmetry.status == MOI.OPTIMAL
@assert isapprox(wedderburn_dec.objective, no_symmetry.objective, atol = 1e-3)
@assert no_symmetry.solve_t / wedderburn_dec.solve_t > 1

