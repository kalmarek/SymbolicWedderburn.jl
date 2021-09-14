function _row_echelon(A; L = [])
    A, l = Characters.row_echelon_form(A)
    @test all([x ∈ l for x in L]) && length(L) <= length(l)
    for (i, j) in enumerate(l)
        @test all([iszero(A[k, j]) for k = 1:size(A, 1) if k != i])
    end
    @test all([isone(A[i, j]) for (i, j) in enumerate(l)])
end

function _nullspace(A)
    V = Characters.left_nullspace(A)
    W = Characters.right_nullspace(A)
    @test all(iszero.(V*A))
    @test all(iszero.(A*W))
end

function _left_eigen(A)
    E = Characters.left_eigen(A)
    for (val, space) in E
        @test space * A == space .* val
    end
end


@testset "RowEchelon" begin
    @testset "examples" begin
        A = GF{11}.([ 0 1 2; 1 2 1; 2 7 8])
        L = [1, 2]
        _row_echelon(A; L = L)

        A = GF{7}.([0 1; 1 2; 0 5])
        L = [1, 2]
        _row_echelon(A; L = L)

        A = GF{2}.([ 0 1 0 0; 0 0 1 1; 0 1 0 1; 0 0 0 0])
        L = [2, 3, 4]
        _row_echelon(A, L = L)

        #Handbook of Computational Group Theory p.252
        A = GF{3}.([1 0 2 0 0 0 0 0;
                    2 1 0 0 0 1 2 0;
                    0 2 1 0 0 0 1 2;
                    0 0 0 0 0 2 0 1])

        N = GF{3}.([1 0 2 0 0 0 0 0;
                    0 1 2 0 0 0 2 1;
                    0 0 0 0 0 1 0 2;
                    0 0 0 0 0 0 0 0])

     @test Characters.row_echelon_form(A)[1] == N

    end
    @testset "random" begin
        for p in [2, 5, 7]
            for s in [3, 5, 6]
                for t in [3, 5, 9]
                    _row_echelon(GF{p}.(rand(Int8, s, t)))
                end
            end
        end
    end
end

@testset "NullSpace" begin
    @testset "examples" begin
        A = [2 1 0 2; 0 0 1 3; 0 1 0 3; 0 0 0 0]
        for p in [2, 3, 5]
            _nullspace(GF{p}.(A))
        end

        A  = [1 2 0 1 4; 0 1 2 1 3; 1 2 1 3 5; 0 1 2 2 7]
        for p in [2, 3, 5]
            _nullspace(GF{p}.(A))
        end

        A = GF{2}.([ 0 1 0 0; 0 0 1 1; 0 1 0 1; 0 0 0 0])
        _nullspace(A)

        #Handbook of Computational Group Theory p.252
        A = GF{3}.([1 1 1 0 0 0 0 0 2 0 2 0;
                    1 1 1 0 0 0 0 0 0 2 0 2;
                    1 1 1 0 0 0 0 0 2 0 2 0;
                    0 0 0 0 0 0 0 0 0 2 0 2;
                    0 0 0 0 0 0 0 0 1 0 1 0;
                    0 0 0 0 0 1 1 1 0 1 0 1;
                    0 0 0 0 0 1 1 1 1 0 1 0;
                    0 0 0 0 0 1 1 1 0 1 0 1])
        N = GF{3}.([1 0 2 0 0 0 0 0;
                    0 1 2 0 0 0 2 1;
                    0 0 0 1 1 0 2 1;
                    0 0 0 0 0 1 0 2])

        K = Characters.left_nullspace(A)
        @test Characters.row_echelon_form(K)[1] == N

   end
    @testset "random" begin
        for p in [2, 5, 7]
            for s in [3, 5, 6]
                for t in [3, 5, 9]
                    _nullspace(GF{p}.(rand(Int8, s, t)))
                end
            end
        end
    end
end

@testset "LeftEigen" begin
    @testset "examples" begin
        #Handbook of Computational Group Theory p.259
        M2 = GF{7}.([0 1 0 0;
                     3 2 0 0;
                     0 0 3 0;
                     0 0 0 3])
        M3 = GF{7}.([0 0 1 0;
                     0 0 3 0;
                     0 0 0 4;
                     4 4 0 0])
        M4 = GF{7}.([0 0 0 1;
                     0 0 0 3;
                     4 4 0 0;
                     0 0 4 0])

        E2 = Characters.left_eigen(M2)
        @test size(E2[GF{7}(3)], 1) == 3

        E3 = Characters.left_eigen(M3)
        E4 = Characters.left_eigen(M4)
        for val in GF{7}.([0, -3, -5, -6])
            @test haskey(E3, val)
            @test haskey(E4, val)
        end

        @test E3[GF{7}(0)] == GF{7}.([1 2 0 0])
        @test E3[GF{7}(-3)] == GF{7}.([1 1 1 1])
        @test E3[GF{7}(-5)] == GF{7}.([1 1 2 4])
        @test E3[GF{7}(-6)] == GF{7}.([1 1 4 2])

    end

    @testset "random" begin
        for p in [2, 5, 7]
            for s in [3, 5, 6, 9]
                _left_eigen(GF{p}.(rand(Int8, s, s)))
            end
        end
    end

end

@testset "EigSpaceDec" begin
      @testset "examples" begin
        #Handbook of Computational Group Theory p.259
        M3 = GF{7}.([0 0 1 0;
                     0 0 3 0;
                     0 0 0 4;
                     4 4 0 0])
        esd = Characters.EigenSpaceDecomposition(M3)
        @test Characters.isdiag(esd)

        #Handbook of Computational Group Theory p.262
        M6 = GF{13}.([0 0 0 0 0 1;
                      0 0 1 0 0 0;
                      0 2 0 0 0 1;
                      0 0 0 0 2 0;
                      0 0 0 2 0 0;
                      2 0 1 0 0 0])
        esd = Characters.EigenSpaceDecomposition(M6)

        @test length(esd) == 4
        @test esd[4] == GF{13}.([1 1 6 0 0 6])
        @test esd[3] == GF{13}.([1 12 1 0 0 12; 0 0 0 1 12 0])
        @test esd[2] == GF{13}.([1 1 1 0 0 1; 0 0 0 1 1 0])
        @test esd[1] == GF{13}.([1 12 6 0 0 7])

        M4 = GF{13}.([0 0 0 1 0 0;
                      0 0 0 0 1 0;
                      0 0 0 2 0 0;
                      0 3 0 0 0 3;
                      3 0 3 0 0 0;
                      0 0 0 0 2 0])

        esd = Characters.refine(esd, M4)
        @test Characters.isdiag(esd)
        @test esd[1] == GF{13}.([1 12 6 0 0 7])
        @test esd[2] == GF{13}.([1 1 1 1 1 1])
        @test esd[3] == GF{13}.([1 1 1 12 12 1])
        @test esd[4] == GF{13}.([1 12 1 8 5 12])
        @test esd[5] == GF{13}.([1 12 1 5 8 12])
        @test esd[6] == GF{13}.([1 1 6 0 0 6])

    end
end
