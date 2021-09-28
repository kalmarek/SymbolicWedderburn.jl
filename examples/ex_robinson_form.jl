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
using SymbolicWedderburn.StarAlgebras
using DynamicPolynomials
using SumOfSquares
using SCS

const OPTIMIZER = optimizer_with_attributes(
    SCS.Optimizer,
    "acceleration_lookback" => 10,
    "max_iters" => 3_000,
    "alpha" => 1.95,
    "eps" => 1e-6,
    "linear_solver" => SCS.DirectSolver,
)

@polyvar x y
const robinson_form = x^6 + y^6 - x^4 * y^2 - y^4 * x^2 - x^4 - y^4 - x^2 - y^2 + 3x^2 * y^2 + 1

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
# definitions of general actions on polynomials:
include(joinpath(@__DIR__, "action_polynomials.jl"))

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

G = DihedralGroup(4)
for g in G
    @assert SymbolicWedderburn.action(DihedralAction(), g, robinson_form) == robinson_form
end

# A simple Gramm-matrix like formulation of SOS problem:
# * a single large pdd constraint and
# * one linear constraint for each monomial in the monomial basis of f.
no_symmetry = let f = robinson_form
    stats = @timed let m = Model(OPTIMIZER)
        @variable m t
        @objective m Max t
        @constraint m f - t in SOSCone()
        m
    end
    m = stats.value
    optimize!(m)
    @info (m,) termination_status(m) objective_value(m) creation_time = (stats.time) solve_time(m)

    @assert isapprox(value(m[:t]), -3825 / 4096, rtol = 1e-2)

    (
        status = termination_status(m),
        objective = objective_value(m),
        creation_t = stats.time,
        solve_t = solve_time(m)
    )
end

# Let's decompose basis into semi-simple summands and prepare the approximately
# reduced optimization problem. Since this is an isomorphism we reduce the
# number of variables from N² to n₁² + … + nₖ², where Σᵢ nᵢ = N. That is instead
# of one psd constraint of size N we have k constraints of sizes n₁,… nₖ.
# The number of linear constraints remains the same as before.
semisimple_dec = let f = robinson_form, G = DihedralGroup(4)
    basis = monomials([x, y], 0:(DynamicPolynomials.maxdegree(f) ÷ 2))

    stats = @timed begin
        sa_basis, symmetry_adaptation_time, = @timed let
            ssimple_basis = SymbolicWedderburn.symmetry_adapted_basis(
                Rational{Int},
                G,
                DihedralAction(),
                basis,
                semisimple = true
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

    @assert isapprox(value(m[:t]), -3825 / 4096, rtol = 1e-2)

    @info (m,) termination_status(m) objective_value(m) solve_time(m) creation_time = (stats.time) symmetry_adaptation_time
    (
        status = termination_status(m),
        objective = objective_value(m),
        creation_t = stats.time,
        solve_t = solve_time(m)
    )
end

# Now given that f must be constant on orbits of monomials under the action of G
# let's exploit that to reduce the number of cnstraints (one for each orbit).
# Here we keep the single large psd constraint of size N×N.
include(joinpath(@__DIR__, "util.jl"))

orbit_dec = let f = robinson_form, G = DihedralGroup(4), T = Float64
    basis_monoms = monomials([x,y], 0:DynamicPolynomials.maxdegree(f))

    invariant_vs, symmetry_adaptation_time, = @timed invariant_vectors(G, DihedralAction(), basis_monoms)

    M = let basis_full = StarAlgebras.Basis{UInt32}(basis_monoms), basis_half = monomials([x,y], 0:(DynamicPolynomials.maxdegree(f) ÷ 2))
        [basis_full[x * y] for x in basis_half, y in basis_half]
    end

    stats = @timed let m = JuMP.Model(OPTIMIZER)
        JuMP.@variable m t
        JuMP.@objective m Max t
        P = JuMP.@variable(m, [1:size(M, 1), 1:size(M, 1)], PSD)

        # preallocating
        M_orb = similar(M, T)

        C = DynamicPolynomials.coefficients(f - t, basis_monoms)

        for iv in invariant_vs
            c = dot(C, iv)

            # invariant_constraint! is defined locally in util.jl
            # it essentially amounts to averaging constraint matrices defined by M
            # with weights provided by iv and storing the result in (dense) matrix M_orb
            M_orb = invariant_constraint!(M_orb, M, iv)

            #=
            if !iszero(M_orb)
                monoms = basis_monoms[unique!(M[map(!iszero, M_orb)])]
                @info "Coefficient of the $monoms: $([DynamicPolynomials.coefficient(f-t, m) for m in monoms])" dot_C_v=c
                @show sparse(M_orb)
            else
                monoms = basis_monoms[unique!(M[SparseArrays.nonzeroinds(iv)])]
                @info "Coefficient of the $monoms: $([DynamicPolynomials.coefficient(f-t, m) for m in monoms])" dot_C_v=c
                @info "M_orb is zero!"
            end
            =#

            JuMP.@constraint m dot(M_orb, P) == c
        end
        m
    end

    m = stats.value
    optimize!(m)

    @assert isapprox(value(m[:t]), -3825 / 4096, rtol = 1e-2)

    @info (m,) termination_status(m) objective_value(m) solve_time(m) creation_time = (stats.time) symmetry_adaptation_time
    (
        status = termination_status(m),
        objective = objective_value(m),
        creation_t = stats.time,
        solve_t = solve_time(m)
    )
end


wedderburn_dec = let f = robinson_form, G = DihedralGroup(4), T = Float64

    basis_half = monomials([x,y], 0:(DynamicPolynomials.maxdegree(f) ÷ 2))

    wedderburn, symmetry_adaptation_time, = @timed WedderburnDecomposition(
        T,
        G,
        DihedralAction(),
        monomials([x,y], 0:DynamicPolynomials.maxdegree(f)),
        basis_half
    )

    M = [SymbolicWedderburn.basis(wedderburn)[x * y] for x in basis_half, y in basis_half]

    stats = @timed let m = JuMP.Model(OPTIMIZER)
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

        C = DynamicPolynomials.coefficients(f - t, SymbolicWedderburn.basis(wedderburn))

        for iv in invariant_vectors(wedderburn)
            c = dot(C, iv)

            # invariant_constraint! is defined locally in util.jl
            # it essentially amounts to averaging constraint matrices defined by M
            # with weights provided by iv and storing the result in (dense) matrix M_orb
            M_orb = invariant_constraint!(M_orb, M, iv)

            #=
            if !iszero(M_orb)
                monoms = SymbolicWedderburn.basis(wedderburn)[unique!(M[map(!iszero, M_orb)])]
                @info "Coefficient of the $monoms: $([DynamicPolynomials.coefficient(f-t, m) for m in monoms])" dot_C_v=c
                @show sparse(M_orb)
            else
                monoms = SymbolicWedderburn.basis(wedderburn)[unique!(M[SparseArrays.nonzeroinds(iv)])]
                @info "Coefficient of the $monoms: $([DynamicPolynomials.coefficient(f-t, m) for m in monoms])" dot_C_v=c
                @info "M_orb is zero!"
            end
            =#

            Mπs = SymbolicWedderburn.diagonalize!(Mπs, M_orb, wedderburn, tmps)

            JuMP.@constraint m sum(
                dot(Mπ, Pπ) for (Mπ, Pπ) in zip(Mπs, psds) if !iszero(Mπ)
            ) == c
        end
        m
    end
    m = stats.value
    optimize!(m)

    @info value(m[:t])
    @assert isapprox(value(m[:t]), -3825 / 4096, rtol = 1e-2)

    @info (m,) termination_status(m) objective_value(m) solve_time(m) creation_time = (stats.time)
    (
        status = termination_status(m),
        objective = objective_value(m),
        creation_t = stats.time,
        solve_t = solve_time(m)
    )
end
