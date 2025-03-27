using Test
using LinearAlgebra, SparseArrays
using SBPOps

@testset "Check mimetic properties for periodic boundary conditions" begin
    N = 64
    h = (1/N)
    grid = 0 : h : 1
    b = ones(N + 1)

    H21, D21, S21_1, S21_N = sbp21d2variable(N + 1, h, b)
    D21 = add_periodic_bcs_SAT21d2(H21, D21, S21_1, S21_N, b)

    v = sin.(2*pi*grid)

    @test ones(N+1)'*H21*v == ones(N+1)'*H21*D21*v
end