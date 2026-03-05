"""Scalar cubic-convolution interpolation matching upstream `cubic_conv_c.c`."""
@inline function _roundval(x::T) where {T<:Real}
    return x > zero(T) ? floor(Int, x + T(0.5)) : -floor(Int, -x + T(0.5))
end

@inline function _k0k1(d::T) where {T<:AbstractFloat}
    if d <= one(T)
        k0 = d * d * (T(2) * d - T(3)) + one(T)
        k1 = d * d * (d - one(T))
        return k0, k1
    elseif d <= T(2)
        k0 = zero(T)
        k1 = d * (d * d - T(5) * d + T(8)) - T(4)
        return k0, k1
    else
        return zero(T), zero(T)
    end
end

@inline function libcconv(a::AbstractMatrix{T}, y::Real, x::Real) where {T}
    ny, nx = size(a)
    F = float(promote_type(typeof(x), typeof(y), real(T)))
    acoef = F(-0.5)

    x_in = F(x)
    y_in = F(y)

    y_round = _roundval(y_in)
    y_pix = clamp(y_round, 2, ny - 2) # 0-based pixel index in C implementation
    yc = F(y_round) - y_in

    ky0_1, ky1_1 = _k0k1(abs(yc - F(2)))
    ky0_2, ky1_2 = _k0k1(abs(yc - F(1)))
    ky0_3, ky1_3 = _k0k1(abs(yc))
    ky0_4, ky1_4 = _k0k1(abs(yc + F(1)))
    ky0_5, ky1_5 = _k0k1(abs(yc + F(2)))

    x_round = _roundval(x_in)
    x_pix = clamp(x_round, 2, nx - 2) # 0-based pixel index in C implementation
    xc = F(x_round) - x_in

    kx0_1, kx1_1 = _k0k1(abs(xc - F(2)))
    kx0_2, kx1_2 = _k0k1(abs(xc - F(1)))
    kx0_3, kx1_3 = _k0k1(abs(xc))
    kx0_4, kx1_4 = _k0k1(abs(xc + F(1)))
    kx0_5, kx1_5 = _k0k1(abs(xc + F(2)))

    acc = zero(promote_type(T, F))

    @inbounds for j in 0:4
        yoff0 = y_pix + j - 2
        (0 <= yoff0 < ny) || continue
        ky0, ky1 = if j == 0
            ky0_1, ky1_1
        elseif j == 1
            ky0_2, ky1_2
        elseif j == 2
            ky0_3, ky1_3
        elseif j == 3
            ky0_4, ky1_4
        else
            ky0_5, ky1_5
        end
        yidx = yoff0 + 1

        for i in 0:4
            xoff0 = x_pix + i - 2
            (0 <= xoff0 < nx) || continue
            kx0, kx1 = if i == 0
                kx0_1, kx1_1
            elseif i == 1
                kx0_2, kx1_2
            elseif i == 2
                kx0_3, kx1_3
            elseif i == 3
                kx0_4, kx1_4
            else
                kx0_5, kx1_5
            end
            k = kx0 * ky0 + acoef * (kx0 * ky1 + kx1 * ky0) + (acoef * acoef) * kx1 * ky1
            acc += k * a[yidx, xoff0 + 1]
        end
    end

    return T <: Real ? real(acc) : acc
end
