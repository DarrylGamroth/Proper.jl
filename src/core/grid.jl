@inline function coordinate_axis(n::Integer, dx::T) where {T<:AbstractFloat}
    c = n ÷ 2
    return (collect(0:(n - 1)) .- c) .* dx
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
