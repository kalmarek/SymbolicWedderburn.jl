using JuMP
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

import COSMO
function cosmo_optimizer(;
    accel = 15,
    alpha = 1.6,
    eps = 1e-8,
    max_iters = 2_000,
    rho = 0.1,
    verbose = true,
    verbose_timing = verbose,
)
    return JuMP.optimizer_with_attributes(
        COSMO.Optimizer,
        "accelerator" => COSMO.with_options(
            COSMO.AndersonAccelerator,
            mem = accel
        ),
        "alpha" => alpha,
        "decompose" => true,
        "eps_abs" => eps,
        "eps_rel" => eps,
        "eps_prim_inf" => eps,
        "eps_dual_inf" => eps,
        "max_iter" => max_iters,
        "rho" => rho,
        "verbose" => verbose,
        "verbose_timing" => verbose_timing,
    )
end
