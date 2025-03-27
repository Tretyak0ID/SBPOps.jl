export norm21, norm42, norm63
export norm21_uw, norm42_uw, norm63_uw

function norm21(N::Int, h::Union{Float32, Float64, Int32, Int64})
    H_vec = ones(N)
    H_vec[1] = 1/2
    H_vec[N] = 1/2
    H_vec = H_vec * h
    H  = spdiagm(0 => H_vec)
    Hi = spdiagm(0 => 1 ./ H_vec)

    return H, Hi
end

function norm42(N::Int, h::Union{Float32, Float64, Int32, Int64})
    H_vec = ones(N)
    H_vec[1:4]   = [17/48, 59/48, 43/48, 49/48]
    H_vec[N-3:N] = [49/48, 43/48, 59/48, 17/48]
    H_vec = H_vec * h
    H  = spdiagm(0 => H_vec)
    Hi = spdiagm(0 => 1 ./ H_vec)

    return H, Hi
end

function norm63(N::Int, h::Union{Float32, Float64, Int32, Int64})
    H_vec = ones(N)
    H_vec[1:6]   = [13649/43200, 12013/8640, 2711/4320, 5359/4320, 7877/8640, 43801/43200]
    H_vec[N-5:N] = [43801/43200, 7877/8640, 5359/4320, 2711/4320, 12013/8640, 13649/43200]
    H_vec = H_vec * h
    H  = spdiagm(0 => H_vec)
    Hi = spdiagm(0 => 1 ./ H_vec)

    return H, Hi
end

function norm21_uw(N::Int, h::Union{Float32, Float64, Int32, Int64})
    H_vec = ones(N)
    H_vec[1:2]   = [1/4, 5/4]
    H_vec[N-1:N] = [5/4, 1/4]
    H_vec = H_vec * h
    H  = spdiagm(0 => H_vec)
    Hi = spdiagm(0 => 1 ./ H_vec)

    return H, Hi
end

function norm42_uw(N::Int, h::Union{Float32, Float64, Int32, Int64})
    H_vec = ones(N)
    H_vec[1:4]   = [49/144, 61/48, 41/48, 149/144]
    H_vec[N-3:N] = [149/144, 41/48, 61/48, 49/144]
    H_vec = H_vec * h
    H  = spdiagm(0 => H_vec)
    Hi = spdiagm(0 => 1 ./ H_vec)

    return H, Hi
end

function norm63_uw(N::Int, h::Union{Float32, Float64, Int32, Int64})
    H_vec = ones(N)
    H_vec[1:6]   = [13613/43200, 12049/8640, 535/864, 1079/864, 7841/8640, 43837/43200]
    H_vec[N-5:N] = [43837/43200, 7841/8640, 1079/864, 535/864, 12049/8640, 13613/43200]
    H_vec = H_vec * h
    H  = spdiagm(0 => H_vec)
    Hi = spdiagm(0 => 1 ./ H_vec)

    return H, Hi
end