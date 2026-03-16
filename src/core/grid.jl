@inline function coordinate_axis(n::Integer, dx::T) where {T<:AbstractFloat}
    return range(-(n ÷ 2) * dx, step=dx, length=n)
end

function radius_map(ny::Integer, nx::Integer, dx::T) where {T<:AbstractFloat}
    y = coordinate_axis(ny, dx)
    x = coordinate_axis(nx, dx)
    r2 = Array{T}(undef, ny, nx)
    @inbounds for j in 1:nx
        x2 = x[j] * x[j]
        for i in 1:ny
            r2[i, j] = sqrt(x2 + y[i] * y[i])
        end
    end
    return r2
end

function spatial_frequency_axis(n::Integer, dx::T) where {T<:AbstractFloat}
    step = inv(T(n) * dx)
    return range(-(n ÷ 2) * step, step=step, length=n)
end

function fft_order_rho2_map(ny::Integer, nx::Integer, dx::T) where {T<:AbstractFloat}
    fy = spatial_frequency_axis(ny, dx)
    fx = spatial_frequency_axis(nx, dx)
    centered = Array{T}(undef, ny, nx)
    @inbounds for j in 1:nx
        fx2 = fx[j] * fx[j]
        for i in 1:ny
            centered[i, j] = fx2 + fy[i] * fy[i]
        end
    end
    return FFTW.ifftshift(centered)
end

function fft_order_rsqr_map(ny::Integer, nx::Integer, dx::T) where {T<:AbstractFloat}
    y = coordinate_axis(ny, dx)
    x = coordinate_axis(nx, dx)
    centered = Array{T}(undef, ny, nx)
    @inbounds for j in 1:nx
        x2 = x[j] * x[j]
        for i in 1:ny
            centered[i, j] = x2 + y[i] * y[i]
        end
    end
    return FFTW.fftshift(centered)
end
