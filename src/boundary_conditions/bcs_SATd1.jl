export add_periodic_SAT_bcs

function add_periodic_SAT_bcs(H::SparseMatrixCSC{<:Real}, D::SparseMatrixCSC{<:Real})
    B = spzeros(size(D));
    B[1,1] = 1/2/H[1,1]; B[1,end] = -1/2/H[1,1];
    B[end,1] = 1/2/H[end,end]; B[end,end] = -1/2/H[end,end];
    D + B
end