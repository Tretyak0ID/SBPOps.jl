using Test
using LinearAlgebra, SparseArrays
using SBPOps

@testset "Check approximation orders test" begin
    for N in (64, 65)
        h = (1/N)
        grid = 0 : h : 1
        b = ones(N + 1)

        H21, D21_2, S21_1, S21_N = sbp21d2(N + 1, h)
        H42, D42_2, S21_2, S21_N = sbp42d2(N + 1, h)
        H21b, D21_2b, S21_1b, S21_Nb = sbp21d2variable(N + 1, h, b)

        @test D21_2 * ones(N + 1) ≈ zeros(N + 1) atol=1e-14
        @test D21_2 * grid ≈ zeros(N + 1) atol=1e-14

        @test D21_2b * ones(N + 1) ≈ zeros(N + 1) atol=1e-14
        @test D21_2b * grid ≈ zeros(N + 1) atol=1e-14

        @test D42_2 * ones(N + 1) ≈ zeros(N + 1) atol=1e-11
        @test D42_2 * grid ≈ zeros(N + 1) atol=1e-11
        @test D42_2 * grid .^ 2 ≈ 2 .* ones(N + 1) atol=1e-10
    end
end

@testset "Check SBP-propertys test" begin
    for N in (64, 65)
        h = (1/N)
        grid = 0 : h : 1
        b = ones(N + 1)

        e1 = spzeros(N + 1, 1); e1[1] = 1
        eN = spzeros(N + 1, 1); eN[N + 1] = 1

        H21, D21 = sbp21(N + 1, h)
        H42, D42 = sbp42(N + 1, h)
        H21, D21_2, S21_1, S21_N = sbp21d2(N + 1, h)
        H42, D42_2, S21_2, S21_N = sbp42d2(N + 1, h)
        H21b, D21_2b, S21_1b, S21_Nb = sbp21d2variable(N + 1, h, b)
        

        u = sin.(grid) + cos.(grid)
        v = 12 * grid.^2

        R21b = Matrix(H21b*D21_2b + D21'*H21b*spdiagm(0=>b)*D21 +  b[1]*e1*S21_1b - b[N]*eN*S21_Nb)
        @test issymmetric(R21b)
        @test minimum(eigvals(R21b)) <= 0
    end
end