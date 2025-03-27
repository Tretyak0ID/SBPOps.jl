using Test

# @testset "sbp_ops.jl" begin
#     # Write your tests here.
# end

@time @testset "SBPOps.jl tests" begin
    @time @testset "Norm for sbp operators test" begin include("norms_tests.jl") end
    @time @testset "Standard sbp operators test" begin include("sbp_standard_tests.jl") end
    @time @testset "Diff2 sbp operators test" begin include("sbp_diff2_tests.jl") end
    @time @testset "SAT boundary conditions for d2-sbp-operators test" begin include("bcs_SAT2d_tests.jl") end
    #@time @testset "Upwind sbp operators test" begin include("sbp_upwind_tests.jl") end
end