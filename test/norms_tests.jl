using Test
using LinearAlgebra, SparseArrays
using SBPOps

@testset "Property of Standard norm operators test" begin
for T in (Float32, Float64), N in (64, 65)
    H21, Hi21 = norm21(N, 1/N)
    H42, Hi42 = norm42(N, 1/N)
    H63, Hi63 = norm63(N, 1/N)

    @test issymmetric(H21) == true
    @test issymmetric(H42) == true
    @test issymmetric(H63) == true
    @test issymmetric(Hi21) == true
    @test issymmetric(Hi42) == true
    @test issymmetric(Hi63) == true
    @test inv(Matrix(H21)) ≈ Hi21
    @test inv(Matrix(H42)) ≈ Hi42
    @test inv(Matrix(H63)) ≈ Hi63

    @test minimum(eigvals(Matrix(H21))) >= 0
    @test minimum(eigvals(Matrix(H42))) >= 0
    @test minimum(eigvals(Matrix(H63))) >= 0
end
end

@testset "Property of Upwind norm operators test" begin
    for T in (Float32, Float64), N in (64, 65)
        H21uw, Hi21uw = norm21_uw(N, 1/N)
        H42uw, Hi42uw = norm42_uw(N, 1/N)
        H63uw, Hi63uw = norm63_uw(N, 1/N)
    
        @test issymmetric(H21uw) == true
        @test issymmetric(H42uw) == true
        @test issymmetric(H63uw) == true
        @test issymmetric(Hi21uw) == true
        @test issymmetric(Hi42uw) == true
        @test issymmetric(Hi63uw) == true
        @test inv(Matrix(H21uw)) ≈ Hi21uw
        @test inv(Matrix(H42uw)) ≈ Hi42uw
        @test inv(Matrix(H63uw)) ≈ Hi63uw
    
        @test minimum(eigvals(Matrix(H21uw))) >= 0
        @test minimum(eigvals(Matrix(H42uw))) >= 0
        @test minimum(eigvals(Matrix(H63uw))) >= 0
    end
    end