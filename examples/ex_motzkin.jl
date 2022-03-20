using SymbolicWedderburn
using PermutationGroups

using DynamicPolynomials

include(joinpath(@__DIR__, "action_polynomials.jl"))
include(joinpath(@__DIR__, "sos_problem.jl"))
include(joinpath(@__DIR__, "solver.jl"))

@polyvar x y

motzkin = x^4 * y^2 + y^4 * x^2 - 3 * x^2 * y^2 + 1
g = (x^2 + y^2 + 1)

no_symmetry = let poly = motzkin * g
    m, model_creation_t = @timed sos_problem(poly)
    JuMP.set_optimizer(
        m,
        scs_optimizer(max_iters = 20_000, eps = 1e-7, accel = 20),
    )

    optimize!(m)
    (
        status = termination_status(m),
        objective = objective_value(m),
        symmetry_adaptation_t = 0.0,
        creation_t = model_creation_t,
        solve_t = solve_time(m),
    )
end

wedderburn_dec =
    let poly = motzkin * g,
        G = PermGroup(perm"(1,2)"),
        action = VariablePermutation()

        m, stats = sos_problem(poly, G, action)
        JuMP.set_optimizer(
            m,
            scs_optimizer(max_iters = 20_000, eps = 1e-7, accel = 20),
        )

        optimize!(m)
        (
            status = termination_status(m),
            objective = objective_value(m),
            symmetry_adaptation_t = stats["symmetry_adaptation"],
            creation_t = stats["model_creation"],
            solve_t = solve_time(m),
        )
    end

@assert no_symmetry.status == wedderburn_dec.status == MOI.OPTIMAL
@assert isapprox(no_symmetry.objective, wedderburn_dec.objective, atol = 1e-6)
@assert no_symmetry.solve_t > wedderburn_dec.solve_t

@info "Summary for motzkin example:" no_symmetry wedderburn_dec
