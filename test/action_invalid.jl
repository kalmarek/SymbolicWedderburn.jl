using SymbolicWedderburn
using DynamicPolynomials

include(joinpath(dirname(@__DIR__), "examples", "action_polynomials.jl"))

@testset "Catching Invalid actions" begin
    include(joinpath(dirname(pathof(GroupsCore)), "..", "test", "cyclic.jl"))
    struct CyclicAction <: OnMonomials end
    SymbolicWedderburn.coeff_type(::CyclicAction) = Float64

    function SymbolicWedderburn.action(
        ::CyclicAction,
        el::CyclicGroupElement,
        mono::AbstractMonomial,
    )
        b = a[2], a[1]
        sign_aodd = 1 <= el.residual <= 2 ? -1 : 1
        sign_aeven = 2 <= el.residual ? -1 : 1
        return mono([a[1], a[2]] => [sign_aodd * b[1], sign_aeven * b[2]])
    end

    G = CyclicGroup(4)
    act = CyclicAction()
    @polyvar a[1:2]
    basis = monomials(a, 0:2)

    @test_throws SymbolicWedderburn.GroupActionError SymbolicWedderburn.check_group_action(
        G,
        act,
        basis,
    )

    function SymbolicWedderburn.action(
        ::CyclicAction,
        el::CyclicGroupElement,
        mono::AbstractMonomial,
    )
        b = el.residual < 2 ? (a[1], a[2]) : (a[2], a[1])
        sign_aodd = 1 <= el.residual <= 2 ? -1 : 1
        sign_aeven = 2 <= el.residual ? -1 : 1
        return mono([a[1], a[2]] => [sign_aodd * b[1], sign_aeven * b[2]])
    end

    @test_throws SymbolicWedderburn.GroupActionError SymbolicWedderburn.check_group_action(
        G,
        act,
        basis,
    )

    @test_throws SymbolicWedderburn.GroupActionError SymbolicWedderburn.symmetry_adapted_basis(
        G,
        act,
        basis,
    )

    @test_throws SymbolicWedderburn.GroupActionError SymbolicWedderburn.WedderburnDecomposition(
        Float64,
        G,
        act,
        basis,
        basis,
    )
end
