export sbp21d2, sbp42d2#, sbp63_2diff
export sbp21d2variable

function sbp21d2(N::Int, h::Union{Float32, Float64, Int32, Int64})
    H, Hi = norm21(N, h)

    D2 = spdiagm(-1 => 1*ones(N-1), 0 => -2*ones(N), 1 => 1*ones(N-1))
    D2[1, 1] =  1
    D2[1, 2] = -2
    D2[1, 3] =  1
    D2[N, N-2] =  1
    D2[N, N-1] = -2
    D2[N, N]   =  1
    D2 = D2 ./ h^2

    S1 = spzeros(1, N); S1[1 : 3]  = [-3/2, 2, -1/2] ./ h
    SN = spzeros(1, N); SN[N-2 : N] = [1/2, -2, 3/2] ./ h
    
    return H, D2, S1, SN
end

function sbp42d2(N::Int, h::Union{Float32, Float64, Int32, Int64})
    H, Hi = norm42(N, h)

    e1 = spzeros(N, 1); e1[1] = 1
    eN = spzeros(N, 1); eN[N] = 1

    M = -spdiagm(-2 => -1/12*ones(N-2), -1 => 16/12*ones(N-1), 0 => -30/12*ones(N), 1 => 16/12*ones(N-1), 2 => -1/12*ones(N-2))
    M_coefs = [9/8 -59/48 1/12 1/48;
              -59/48 59/24 -59/48 0;
              1/12 -59/48 55/24 -59/48;
              1/48 0 -59/48 59/24]
    M[1:4, 1:4] = M_coefs
    M[N-3:N, N-3:N] = rotr90(M_coefs, 2)
    M = M ./ h

    S1 = spzeros(1, N); S1[1 : 4]   = [-11/6, 3, -3/2, 1/3] ./ h
    SN = spzeros(1, N); SN[N-3 : N] = [-1/3, 3/2, -3, 11/6] ./ h

    D2 = Hi * (-M - e1*S1 + eN*SN)
    
    return H, D2, S1, SN
end

#diff2 sbp variables operators

function sbp21d2variable(N::Int, h::Union{Float32, Float64, Int32, Int64}, b::AbstractVector)

    if (length(b) != N)
        println("ERROR! Dimension veriable coefficients vector is not equal to N!")
        return 0
    end

    H, Hi = norm21(N, h)
    e1 = spzeros(N, 1); e1[1] = 1
    eN = spzeros(N, 1); eN[N] = 1
    S1 = spzeros(1, N); S1[1 : 3]  = [-3/2, 2, -1/2] ./ h
    SN = spzeros(1, N); SN[N-2 : N] = [1/2, -2, 3/2] ./ h

    M_coefs_lu = [(b[1]+b[2])/2 -(b[1]+b[2])/2 0;
                -(b[1]+b[2])/2 b[1]/2+b[2]+b[3]/2 -(b[2]+b[3])/2;
                0 -(b[2]+b[3])/2 b[2]/2+b[3]+b[4]/2]
    M_coefs_rd = [b[N-1]/2+b[N-2]+b[N-3]/2 -(b[N-1]+b[N-2])/2 0;
                -(b[N-1]+b[N-2])/2 b[N]/2+b[N-1]+b[N-2]/2 -(b[N]+b[N-1])/2;
                0 -(b[N]+b[N-1])/2 (b[N]+b[N-1])/2]
    ml = -(b[1:N-1] + b[2:N]) / 2
    mc = ones(N); mc[2:N-1] = (b[1:N-2] + 2*b[2:N-1] + b[3:N]) / 2
    mr = -(b[1:N-1] + b[2:N]) / 2

    M = spdiagm(-1 => ml,  0 => mc, 1 => mr)
    M[1:3, 1:3]     = M_coefs_lu
    M[N-2:N, N-2:N] = M_coefs_rd
    M = M ./ h

    D2 = Hi * (-M - b[1]*e1*S1 + b[N]*eN*SN)

    return H, D2, S1, SN
end