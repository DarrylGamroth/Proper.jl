"""Scalar cubic-convolution interpolation with clamped boundaries."""
@inline function _cubic_weight(t::T) where {T<:AbstractFloat}
    a = T(-0.5)
    u = abs(t)
    if u < one(T)
        return (a + 2) * u^3 - (a + 3) * u^2 + one(T)
    elseif u < T(2)
        return a * u^3 - 5a * u^2 + 8a * u - 4a
    else
        return zero(T)
    end
end

@inline function libcconv(a::AbstractMatrix{T}, y::Real, x::Real) where {T}
    ny, nx = size(a)
    xf = float(x)
    yf = float(y)

    x0 = floor(Int, xf)
    y0 = floor(Int, yf)

    acc = zero(T)
    @inbounds for m in -1:2
        yy = _clamp_index(y0 + m, ny)
        wy = _cubic_weight(yf - (y0 + m))
        for n in -1:2
            xx = _clamp_index(x0 + n, nx)
            wx = _cubic_weight(xf - (x0 + n))
            acc += (wy * wx) * a[yy, xx]
        end
    end
    return acc
end
