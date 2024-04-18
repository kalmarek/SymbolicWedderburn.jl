using LinearAlgebra
using SparseArrays

using SymbolicWedderburn
using DynamicPolynomials

include(joinpath(dirname(@__DIR__), "examples", "action_polynomials.jl"))
include(joinpath(dirname(@__DIR__), "examples", "dihedral.jl"))
include(joinpath(dirname(@__DIR__), "examples", "sos_problem.jl"))

using JuMP
import SCS

function scs_optimizer(;
    accel = 0,
    alpha = 1.5,
    eps = 1e-6,
    max_iters = 10_000,
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
        "warm_start" => true,
        "verbose" => verbose,
    )
end

@polyvar x y
const robinson_form =
    x^6 + y^6 - x^4 * y^2 - y^4 * x^2 - x^4 - y^4 - x^2 - y^2 + 3x^2 * y^2 + 1

struct DihedralAction <: OnMonomials end

SymbolicWedderburn.coeff_type(::DihedralAction) = Float64
function SymbolicWedderburn.action(
    ::DihedralAction,
    el::DihedralElement,
    mono::AbstractMonomial,
)
    if iseven(el.reflection + el.id)
        var_x, var_y = x, y
    else
        var_x, var_y = y, x
    end
    sign_x = 1 <= el.id <= 2 ? -1 : 1
    sign_y = 2 <= el.id ? -1 : 1
    return mono([x, y] => [sign_x * var_x, sign_y * var_y])
end

function StarAlgebras.comparable(::Type{DihedralElement})
    return StarAlgebras.Comparable((a, b) -> hash(a) < hash(b))
end

@testset "Dihedral Action" begin
    G = DihedralGroup(4)
    @test all(
        SymbolicWedderburn.action(DihedralAction(), g, robinson_form) ==
        robinson_form for g in DihedralGroup(4)
    )

    T = Float64

    m, _ = sos_problem(robinson_form, G, DihedralAction())
    JuMP.set_optimizer(
        m,
        scs_optimizer(; eps = 1e-5, alpha = 1.95, accel = -15),
    )

    optimize!(m)

    @test isapprox(value(m[:t]), -3825 / 4096, rtol = 1e-4)
    status = termination_status(m)
    @test status ∈ (MOI.OPTIMAL, MOI.ALMOST_OPTIMAL)
end

# same action but through BySignedPermutations:
struct DihedralActionSP <: SymbolicWedderburn.BySignedPermutations end

function SymbolicWedderburn.action(
    ::DihedralActionSP,
    el::DihedralElement,
    mono::AbstractMonomial,
)
    var_x, var_y = iseven(el.reflection + el.id) ? (x, y) : (y, x)

    sign_x = 1 <= el.id <= 2 ? -1 : 1
    sign_y = 2 <= el.id ? -1 : 1

    g_mono = mono([x, y] => [sign_x * var_x, sign_y * var_y])
    sign = DP.coefficient(mono) ÷ DP.coefficient(g_mono)

    return monomial(g_mono), sign
end

# This is only needed to define action on the whole Robinson form to check that is actually invariant.
function SymbolicWedderburn.action(
    ac::SymbolicWedderburn.BySignedPermutations,
    g::DihedralElement,
    poly::AbstractPolynomial,
)
    return sum(monomials(poly)) do m
        c = DynamicPolynomials.coefficient(poly, m)
        gm, u = SymbolicWedderburn.action(ac, g, m)
        return u * c * gm
    end
end

@testset "Dihedral action through Signed Permutations" begin
    @test all(
        SymbolicWedderburn.action(DihedralAction(), g, mono) ==
        prod(SymbolicWedderburn.action(DihedralActionSP(), g, mono)) for
        mono in monomials([x, y], 0:4), g in DihedralGroup(4)
    )

    m, _ = sos_problem(robinson_form, DihedralGroup(4), DihedralActionSP())

    JuMP.set_optimizer(
        m,
        scs_optimizer(;
            max_iters = 5_000,
            alpha = 1.8,
            accel = -15,
            eps = 1e-5,
        ),
    )
    optimize!(m)
    @test termination_status(m) == MOI.OPTIMAL
    @test isapprox(objective_value(m), -3825 / 4096, atol = 1e-3)

    @test all(
        SymbolicWedderburn.action(DihedralActionSP(), g, robinson_form) ==
        robinson_form for g in DihedralGroup(4)
    )
end
