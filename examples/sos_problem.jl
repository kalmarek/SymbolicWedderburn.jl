using LinearAlgebra
using SparseArrays

using GroupsCore

using DynamicPolynomials
using SymbolicWedderburn
import SymbolicWedderburn.StarAlgebras

using JuMP

function invariant_constraint!(
    M_orb::AbstractMatrix{<:AbstractFloat},
    M::Matrix{<:Integer},
    invariant_vec::SparseVector,
)
    M_orb .= zero(eltype(M_orb))
    for i in eachindex(M)
        if M[i] ∈ SparseArrays.nonzeroinds(invariant_vec)
            M_orb[i] += invariant_vec[M[i]]
        end
    end
    return M_orb
end

function sos_problem(poly::AbstractPolynomial)
    vars = DynamicPolynomials.variables(poly)

    basis_psd = DynamicPolynomials.monomials(
        vars,
        0:DynamicPolynomials.maxdegree(poly)÷2,
    )

    basis_constraints = StarAlgebras.Basis{UInt16}(
        DynamicPolynomials.monomials(
            vars,
            0:DynamicPolynomials.maxdegree(poly),
        ),
    )

    M = [basis_constraints[x*y] for x in basis_psd, y in basis_psd]

    sos_model = JuMP.Model()
    JuMP.@variable sos_model t
    JuMP.@objective sos_model Max t
    n = length(basis_psd)
    P = JuMP.@variable sos_model P[1:n, 1:n] Symmetric
    JuMP.@constraint sos_model P in PSDCone()

    objective = poly - t
    for (idx, b) in enumerate(basis_constraints)
        c = DynamicPolynomials.coefficient(objective, b)
        JuMP.@constraint sos_model LinearAlgebra.dot(P, M .== idx) == c
    end

    return sos_model
end

function sos_problem(
    poly::AbstractPolynomial,
    invariant_vs::AbstractVector,
    basis_constraints::StarAlgebras.AbstractBasis,
    basis_psd,
    T = Float64,
)
    M = [basis_constraints[x*y] for x in basis_psd, y in basis_psd]

    sos_model = JuMP.Model()
    JuMP.@variable sos_model t
    JuMP.@objective sos_model Max t
    n = length(basis_psd)
    P = JuMP.@variable sos_model P[1:n, 1:n] Symmetric
    JuMP.@constraint sos_model P in PSDCone()

    # preallocating
    M_orb = similar(M, T)

    C = DynamicPolynomials.coefficients(poly - t, basis_constraints)

    for iv in invariant_vs
        c = dot(C, iv)
        # average Ms into M_orb with weights given by iv
        M_orb = invariant_constraint!(M_orb, M, iv)
        JuMP.@constraint sos_model dot(M_orb, P) == c
    end

    return sos_model
end

function sos_problem(
    poly::AbstractPolynomial,
    wedderburn::SymbolicWedderburn.WedderburnDecomposition,
    basis_psd;,
)
    m = JuMP.Model()

    M = let basis_constraints = SymbolicWedderburn.basis(wedderburn)
        [basis_constraints[x*y] for x in basis_psd, y in basis_psd]
    end

    JuMP.@variable m t
    JuMP.@objective m Max t
    psds = map(SymbolicWedderburn.direct_summands(wedderburn)) do ds
        dim = size(ds, 1)
        P = JuMP.@variable m [1:dim, 1:dim] Symmetric
        JuMP.@constraint m P in PSDCone()
        return P
    end

    # preallocating
    Mπs = zeros.(eltype(wedderburn), size.(psds))
    M_orb = similar(M, eltype(wedderburn))
    tmps = SymbolicWedderburn._tmps(wedderburn)

    C = DynamicPolynomials.coefficients(
        poly - t,
        SymbolicWedderburn.basis(wedderburn),
    )
    for iv in invariant_vectors(wedderburn)
        c = dot(C, iv)
        M_orb = invariant_constraint!(M_orb, M, iv)
        Mπs = SymbolicWedderburn.diagonalize!(Mπs, M_orb, wedderburn, tmps)

        JuMP.@constraint m sum(
            dot(Mπ, Pπ) for (Mπ, Pπ) in zip(Mπs, psds) if !iszero(Mπ)
        ) == c
    end
    return m
end

function sos_problem(
    poly::AbstractPolynomial,
    G::Group,
    action::SymbolicWedderburn.Action,
    T = Float64;
    decompose_psd = true,
    semisimple = false,
)
    max_deg = DynamicPolynomials.maxdegree(poly)
    vars = DynamicPolynomials.variables(poly)
    basis_psd = DynamicPolynomials.monomials(vars, 0:max_deg÷2)
    basis_constraints = DynamicPolynomials.monomials(vars, 0:max_deg)

    if decompose_psd == true
        wedderburn, symmetry_adaptation_time = @timed WedderburnDecomposition(
            T,
            G,
            action,
            basis_constraints,
            basis_psd,
            semisimple = semisimple,
        )

        model, model_creation_time =
            @timed sos_problem(poly, wedderburn, basis_psd)
    else
        (invariant_vs, basis_cnstr), symmetry_adaptation_time = @timed let G = G
            basis = StarAlgebras.Basis{UInt32}(basis_constraints)

            tblG = SymbolicWedderburn.Characters.CharacterTable(
                Rational{Int},
                G,
            )
            iv = invariant_vectors(tblG, action, basis)
            iv, basis
        end

        model, model_creation_time =
            @timed sos_problem(poly, invariant_vs, basis_cnstr, basis_psd, T)
    end

    stats = Dict(
        "symmetry_adaptation" => symmetry_adaptation_time,
        "model_creation" => model_creation_time,
    )

    return model, stats
end
