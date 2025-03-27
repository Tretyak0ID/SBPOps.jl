export add add_periodic_SAT_bcs

function add_periodic_SAT_bcs(H::SparseMatrixCSC{<:Real}, D::SparseMatrixCSC{<:Real})
    B = zeros(size(D))