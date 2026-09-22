import SCS

function scs_optimizer(;
    accel = 0,
    alpha = 1.5,
    eps = 1e-6,
    max_iters = 100_000,
    rho = 1e-6,
    verbose = true,
)
    # Preserve the pre-3.3 acceleration and scaling behavior.
    settings = if VersionNumber(SCS.scs_version()) >= v"3.3.0"
        [
            "acceleration_lookback" => abs(accel),
            "acceleration_type_1" => Int(accel > 0),
            "acceleration_regularization" => accel > 0 ? 1e-8 : 1e-12,
            "adaptive_diag_scale" => false,
        ]
    else
        ["acceleration_lookback" => accel]
    end
    return JuMP.optimizer_with_attributes(
        SCS.Optimizer,
        settings...,
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
