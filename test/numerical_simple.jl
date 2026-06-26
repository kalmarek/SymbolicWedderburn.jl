# Tests for the numerical block-diagonalization primitives that were ported
# from `SumOfSquares/Certificate/Symmetry/block_diag.jl`.

function _test_orthogonal_transformation_to(A, B, As, Bs)
    U = SW.orthogonal_transformation_to(A, B, As, Bs)
    @test A ≈ U' * B * U
end

function _test_orthogonal_transformation_to(λ, As, Bs)
    A = sum(λ[i] * As[i] for i in eachindex(As))
    B = sum(λ[i] * Bs[i] for i in eachindex(Bs))
    _test_orthogonal_transformation_to(A, B, As, Bs)
    return
end

function _test_orthogonal_transformation_to(A, B)
    _test_orthogonal_transformation_to(A, B, typeof(A)[], typeof(B)[])
    _test_orthogonal_transformation_to(A, B, [A], [B])
    return
end

function _test_orthogonal_transformation_to(T::Type)
    A1 = T[
        0 -1
        1 0
    ]
    A2 = T[
        0 1
        -1 0
    ]
    D = Diagonal(T[1, -1])
    @test A1 == D * A2 * D
    _test_orthogonal_transformation_to(A1, A2)
    B1 = T[
        0 1
        1 0
    ]
    B2 = T[
        0 -1
        -1 0
    ]
    @test B1 == D * B2 * D
    _test_orthogonal_transformation_to(B1, B2)
    @test (A1 + B1) == D * (A2 + B2) * D
    _test_orthogonal_transformation_to(A1 + B1, A2 + B2)
    A1 = T[
        0 1
        -2 0
    ]
    A2 = T[
        0 -1
        2 0
    ]
    _test_orthogonal_transformation_to(A1, A2)
    A1 = T[
        0 -1
        -2 0
    ]
    A2 = T[
        0 1
        2 0
    ]
    _test_orthogonal_transformation_to(A1, A2)
    A1 = T[
        0 0 1
        1 0 0
        0 1 0
    ]
    A2 = T[
        0 1 0
        0 0 1
        1 0 0
    ]
    _test_orthogonal_transformation_to(A1, A2)
    A1 = T[
        1 0 0
        0 -1 0
        0 0 -1
    ]
    A2 = T[
        -1 0 0
        0 -1 0
        0 0 1
    ]
    _test_orthogonal_transformation_to(A1, A2)
    A1 = T[
        -1 0 1
        1 1 0
        0 1 -1
    ]
    A2 = T[
        -1 1 0
        0 1 1
        1 0 -1
    ]
    _test_orthogonal_transformation_to(A1, A2)
    A1 = T[
        -1 1 1
        2 1 0
        0 1 0
    ]
    A2 = T[
        0 1 0
        0 1 2
        1 1 -1
    ]
    _test_orthogonal_transformation_to(A1, A2)
    A1 = T[
        0 1 1
        2 0 0
        0 1 -1
    ]
    A2 = T[
        -1 1 0
        0 0 2
        1 1 0
    ]
    _test_orthogonal_transformation_to(A1, A2)
    A1 = ComplexF64[
        0 1 1
        2 0 0
        0 1 -1
    ]
    A2 = ComplexF64[
        -1 1 0
        0 0 2
        1 1 0
    ]
    _test_orthogonal_transformation_to(A1, A2)
    A1 = T[1 0; 0 -1]
    A2 = T[0 1; 1 0]
    _test_orthogonal_transformation_to([1, 1], [A1, A2], [-A1, A2])
    _test_orthogonal_transformation_to([2, 1], [A1, A2], [-A1, A2])
    return
end

function _test_block_diag(A, d)
    U = SW.ordered_block_diagonalize(A, d)
    @test SW._ordered_block_check(U, A, d)
    return
end

@testset "orthogonal_transformation_to" begin
    @testset "$T" for T in [Int, Float64, ComplexF64]
        _test_orthogonal_transformation_to(T)
    end
end

@testset "ordered_block_diagonalize (dihedral)" begin
    # From SumOfSquares' `dihedral.jl` example
    d = 2
    A2 = [
        0 1 0 0 0 0
        1 0 0 0 0 0
        0 0 0 0 0 1
        0 0 0 0 1 0
        0 0 0 1 0 0
        0 0 1 0 0 0
    ]
    A1 = [
        0 -1 0 0 0 0
        1 0 0 0 0 0
        0 0 0 0 0 -1
        0 0 0 0 1 0
        0 0 0 -1 0 0
        0 0 1 0 0 0
    ]
    _test_block_diag([A1, A2], d)
    # Using `GroupsCore.gens(G::DihedralGroup) = [DihedralElement(G.n, true, 1), DihedralElement(G.n, true, 0)]`, we get
    # see https://github.com/jump-dev/SumOfSquares.jl/issues/381#issuecomment-2296711306
    A1 = Matrix(Diagonal([1, -1, 1, -1, 1, -1]))
    _test_block_diag([A1, A2], d)
end

@testset "ordered_block_diagonalize (alpha)" begin
    α = 0.75
    A = [
        Matrix{Int}(I, 6, 6),
        [
            1 0 0 0 0 0
            0 -1 0 0 0 0
            0 0 0 0 1 0
            0 0 0 0 0 -1
            0 0 1 0 0 0
            0 0 0 -1 0 0
        ],
        [
            -0.5 α 0 0 0 0
            -α -0.5 0 0 0 0
            0 0 -0.5 α 0 0
            0 0 -α -0.5 0 0
            0 0 0 0 -0.5 α
            0 0 0 0 -α -0.5
        ],
    ]
    d = 2
    _test_block_diag(A, d)
end
