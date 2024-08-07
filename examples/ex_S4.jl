const N = 4
@polyvar x[1:N]

OPTIMIZER = scs_optimizer(; max_iters = 2_000, accel = 10, eps = 1e-7)

f =
    1 +
    sum(x .+ 1) +
    sum((x .+ 1) .^ 2)^4 +
    sum((x .+ x') .^ 2)^2 * sum((x .+ 1) .^ 2)

no_symmetry = let f = f
    m, creation_time = @timed sos_problem(f)
    JuMP.set_optimizer(m, OPTIMIZER)
    optimize!(m)

    (
        status = termination_status(m),
        objective = objective_value(m),
        symmetry_adaptation_t = 0.0,
        creation_t = creation_time,
        solve_t = solve_time(m),
    )
end

orbit_dec = let f = f, vars = DynamicPolynomials.variables(f)
    m, stats = sos_problem(
        f,
        PermGroup([perm"(1,2)", Perm([2:N; 1])]),
        VariablePermutation(x);
        decompose_psd = false,
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

semisimple_dec = let f = f, vars = DynamicPolynomials.variables(f)
    mdeg = DynamicPolynomials.maxdegree(f)
    m, stats = sos_problem(
        f,
        PermGroup([perm"(1,2)", Perm([2:N; 1])]),
        VariablePermutation(x);
        semisimple = true,
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

wedderburn_dec = let f = f, vars = DynamicPolynomials.variables(f)
    mdeg = DynamicPolynomials.maxdegree(f)
    m, stats = sos_problem(
        f,
        PermGroup([perm"(1,2)", Perm([2:N; 1])]),
        VariablePermutation(x),
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


@info "Summary of Example: S4" no_symmetry orbit_dec semisimple_dec wedderburn_dec

@assert wedderburn_dec.status == no_symmetry.status == MOI.OPTIMAL
@assert isapprox(wedderburn_dec.objective, no_symmetry.objective, rtol = 1e-5)
@assert no_symmetry.solve_t > wedderburn_dec.solve_t
@assert no_symmetry.solve_t / wedderburn_dec.solve_t > 5
