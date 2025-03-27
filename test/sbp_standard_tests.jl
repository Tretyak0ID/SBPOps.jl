using Test
using LinearAlgebra, SparseArrays
using SBPOps

@testset "Check approximation order test" begin
    for N in (64, 65)
        h = (1/N)
        grid = 0 : h : 1

        H21, D21 = sbp21(N + 1, h)
        H42, D42 = sbp42(N + 1, h)

        @test D21 * ones(N + 1) ≈ zeros(N + 1) atol=1e-14
        @test D21 * grid ≈ ones(N + 1) atol=1e-14

        @test D42 * ones(N + 1) ≈ zeros(N + 1) atol=1e-13
        @test D42 * grid ≈ ones(N + 1) atol=1e-13
        @test D42 * grid .^ 2 ≈ 2 .* grid atol=1e-13
    end
end

@testset "Check SBP-property test" begin
    for N in (64, 65)
        h = (1/N)
        grid = 0 : h : 1

        H21, D21 = sbp21(N + 1, h)
        H42, D42 = sbp42(N + 1, h)
        Q = diagm(zeros(N + 1))
        Q[1, 1]         = -1
        Q[N + 1, N + 1] =  1

        @test Matrix(H21 * D21 + D21' * H21) ≈ Q atol = 1e-13
        @test Matrix(H42 * D42 + D42' * H42) ≈ Q atol = 1e-13
    end
end