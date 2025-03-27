module SBPOps

    using SparseArrays

    include("norms.jl")

    include("operators/sbp_standard.jl")
    include("operators/sbp_upwind.jl")
    include("operators/sbp_diff2.jl")

    include("boundary_conditions/bcs_SATd2.jl")
    include("boundary_conditions/bcs_SATd1.jl")
end
