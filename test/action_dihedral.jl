using LinearAlgebra
using SparseArrays

using SymbolicWedderburn
using DynamicPolynomials

include(joinpath(dirname(@__DIR__), "examples", "action_polynomials.jl"))
include(joinpath(dirname(@__DIR__), "examples", "dihedral.jl"))
include(joinpath(dirname(@__DIR__), "examples", "util.jl"))

using JuMP
using SCS

const OPTIMIZER = optimizer_with_attributes(
    SCS.Optimizer,
    "acceleration_lookback" => 10,
    "max_iters" => 3_000,
    "alpha" => 1.2,
    "eps" => 1e-6,
    "linear_solver" => SCS.DirectSolver,
)

@polyvar x y
const robinson_form = x^6 + y^6 - x^4 * y^2 - y^4 * x^2 - x^4 - y^4 - x^2 - y^2 + 3x^2 * y^2 + 1

struct DihedralAction <: OnMonomials end

SymbolicWedderburn.coeff_type(::DihedralAction) = Float64
function SymbolicWedderburn.action(::DihedralAction, el::DihedralElement, mono::AbstractMonomial)
    if iseven(el.reflection + el.id)
        var_x, var_y = x, y
    else
        var_x, var_y = y, x
    end
    sign_x = 1 <= el.id <= 2 ? -1 : 1
    sign_y = 2 <= el.id ? -1 : 1
    return mono([x, y] => [sign_x * var_x, sign_y * var_y])
end

@testset "Dihedral Action" begin
    G = DihedralGroup(4)
    @test all(
        SymbolicWedderburn.action(DihedralAction(), g, robinson_form) == robinson_form
        for g in DihedralGroup(4)
    )

    T = Float64

    wedderburn = WedderburnDecomposition(
        T,
        G,
        DihedralAction(),
        monomials([x,y], 0:DynamicPolynomials.maxdegree(robinson_form)),
        monomials([x,y], 0:(DynamicPolynomials.maxdegree(robinson_form) ÷ 2))
    )

    M = let bh = monomials([x,y], 0:(DynamicPolynomials.maxdegree(robinson_form) ÷ 2))
        [SymbolicWedderburn.basis(wedderburn)[x * y] for x in bh, y in bh]
    end

    m = let m = JuMP.Model(OPTIMIZER)
        JuMP.@variable m t
        JuMP.@objective m Max t
        psds = [
            JuMP.@variable(m, [1:d, 1:d] in PSDCone())
            for d in size.(direct_summands(wedderburn), 1)
        ]

        # preallocating
        M_orb = similar(M, T)
        # these could be preallocated as well:
        # Mπs = zeros.(T, size.(psds))
        # tmps = SymbolicWedderburn._temps(wedderburn)

        C = DynamicPolynomials.coefficients(robinson_form - t, basis(wedderburn))

        for iv in invariant_vectors(wedderburn)
            c = dot(C, iv)
            M_orb = invariant_constraint!(M_orb, M, iv)
            # Mπs = SymbolicWedderburn.diagonalize!(Mπs, M_orb, wedderburn, tmps)
            Mπs = SymbolicWedderburn.diagonalize(M_orb, wedderburn)

            JuMP.@constraint m sum(
                dot(Mπ, Pπ) for (Mπ, Pπ) in zip(Mπs, psds) if !iszero(Mπ)
            ) == c
        end
        m
    end

    optimize!(m)

    @test isapprox(value(m[:t]), -3825 / 4096, rtol = 1e-4)
    status = termination_status(m)
    @test status ∈ (MOI.OPTIMAL, MOI.ALMOST_OPTIMAL)
end
