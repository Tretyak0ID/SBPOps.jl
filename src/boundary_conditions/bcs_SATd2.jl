export add_periodic_bcs_SAT21d2

function add_periodic_bcs_SAT21d2(H::SparseMatrixCSC{<:Real}, D::SparseMatrixCSC{<:Real}, S1::SparseMatrixCSC{<:Real}, SN::SparseMatrixCSC{<:Real}, b::AbstractVector)
    N = length(S1)
    calc_df_op = spzeros(1, N); calc_df_op[1] = -1; calc_df_op[N] = 1 #calc fN - f1 value for [f1,...,fN]' grid-function 
    extend_op = spzeros(N, 1); extend_op[1] = 1; extend_op[N] = 1     #extend constant alpha of interfaces [alpha,0,...,0,alpha]'
    boundary_diff = spzeros(N,N); boundary_diff[1,:] = -S1; boundary_diff[N,:] = SN
    calc_bound_sum_op = spzeros(1, N); calc_bound_sum_op[1] = 1; calc_bound_sum_op[N] = 1

    SAT_op = 1/2*(spdiagm(0=>b) * boundary_diff)' * extend_op * calc_df_op
    SAT_op = SAT_op - 1/2*extend_op * calc_bound_sum_op * spdiagm(0=>b) * boundary_diff

    return D + SAT_op
end