import SCS

function scs_optimizer(;
    accel = 0,
    alpha = 1.5,
    eps = 1e-6,
    max_iters = 100_000,
    rho = 1e-6,
    verbose = true,
)
    return JuMP.optimizer_with_attributes(
        SCS.Optimizer,
        "acceleration_lookback" => accel,
        "acceleration_interval" => 10,
        "alpha" => alpha,
        "eps_abs" => eps,
        "eps_rel" => eps,
        "linear_solver" => SCS.DirectSolver,
        "max_iters" => max_iters,
        "rho_x" => rho,
        "warm_start" => true,
        "verbose" => verbose,
    )
end
