using KernelAbstractions
import AcceleratedKernels as AK

@inline function _ka_roundval(x::T) where {T<:Real}
    return x > zero(T) ? floor(Int, x + T(0.5)) : -floor(Int, -x + T(0.5))
end

@inline function _ka_k0k1(d::T) where {T<:AbstractFloat}
    if d <= one(T)
        k0 = d * d * (T(2) * d - T(3)) + one(T)
        k1 = d * d * (d - one(T))
        return k0, k1
    elseif d <= T(2)
        k0 = zero(T)
        k1 = d * (d * d - T(5) * d + T(8)) - T(4)
        return k0, k1
    end
    return zero(T), zero(T)
end

@inline function _ka_clamp_index(i::Int, n::Int)
    return ifelse(i < 1, 1, ifelse(i > n, n, i))
end

@inline function _ka_bilinear_sample(a::AbstractMatrix{T}, y::Real, x::Real) where {T}
    ny, nx = size(a)
    x0 = floor(Int, x)
    y0 = floor(Int, y)
    x1 = x0 + 1
    y1 = y0 + 1

    x0c = _ka_clamp_index(x0, nx)
    x1c = _ka_clamp_index(x1, nx)
    y0c = _ka_clamp_index(y0, ny)
    y1c = _ka_clamp_index(y1, ny)

    tx = x - x0
    ty = y - y0

    v00 = a[y0c, x0c]
    v10 = a[y0c, x1c]
    v01 = a[y1c, x0c]
    v11 = a[y1c, x1c]

    v0 = (one(tx) - tx) * v00 + tx * v10
    v1 = (one(tx) - tx) * v01 + tx * v11
    return (one(ty) - ty) * v0 + ty * v1
end

@inline function _ka_cubic_sample(a::AbstractMatrix{T}, y::Real, x::Real) where {T}
    ny, nx = size(a)
    F = float(promote_type(typeof(x), typeof(y), real(T)))
    acoef = F(-0.5)

    x_in = F(x)
    y_in = F(y)

    y_round = _ka_roundval(y_in)
    y_pix = clamp(y_round, 2, ny - 2)
    yc = F(y_round) - y_in

    ky0_1, ky1_1 = _ka_k0k1(abs(yc - F(2)))
    ky0_2, ky1_2 = _ka_k0k1(abs(yc - F(1)))
    ky0_3, ky1_3 = _ka_k0k1(abs(yc))
    ky0_4, ky1_4 = _ka_k0k1(abs(yc + F(1)))
    ky0_5, ky1_5 = _ka_k0k1(abs(yc + F(2)))

    x_round = _ka_roundval(x_in)
    x_pix = clamp(x_round, 2, nx - 2)
    xc = F(x_round) - x_in

    kx0_1, kx1_1 = _ka_k0k1(abs(xc - F(2)))
    kx0_2, kx1_2 = _ka_k0k1(abs(xc - F(1)))
    kx0_3, kx1_3 = _ka_k0k1(abs(xc))
    kx0_4, kx1_4 = _ka_k0k1(abs(xc + F(1)))
    kx0_5, kx1_5 = _ka_k0k1(abs(xc + F(2)))

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

@kernel function _ka_apply_shifted_mask_kernel!(
    field,
    mask,
    sy::Int,
    sx::Int,
    ny::Int,
    nx::Int,
    invert::Bool,
)
    I = @index(Global, NTuple)
    i = I[1]
    j = I[2]

    if i <= ny && j <= nx
        is = i + sy
        js = j + sx
        is = ifelse(is > ny, is - ny, is)
        js = ifelse(js > nx, js - nx, js)
        m = mask[is, js]
        field[i, j] *= invert ? (one(m) - m) : m
    end
end

@kernel function _ka_copy_shifted_complex_kernel!(
    out,
    field,
    r0::Int,
    c0::Int,
    sy::Int,
    sx::Int,
    ny::Int,
    nx::Int,
    oy::Int,
    ox::Int,
)
    I = @index(Global, NTuple)
    i = I[1]
    j = I[2]

    if i <= oy && j <= ox
        is = r0 + i - 1 + sy
        js = c0 + j - 1 + sx
        is = ifelse(is > ny, is - ny, is)
        js = ifelse(js > nx, js - nx, js)
        out[i, j] = field[is, js]
    end
end

@kernel function _ka_copy_shifted_intensity_kernel!(
    out,
    field,
    r0::Int,
    c0::Int,
    sy::Int,
    sx::Int,
    ny::Int,
    nx::Int,
    oy::Int,
    ox::Int,
)
    I = @index(Global, NTuple)
    i = I[1]
    j = I[2]

    if i <= oy && j <= ox
        is = r0 + i - 1 + sy
        js = c0 + j - 1 + sx
        is = ifelse(is > ny, is - ny, is)
        js = ifelse(js > nx, js - nx, js)
        out[i, j] = abs2(field[is, js])
    end
end

@kernel function _ka_cubic_conv_grid_kernel!(
    out,
    a,
    xval,
    yval,
    oy::Int,
    ox::Int,
)
    I = @index(Global, NTuple)
    i = I[1]
    j = I[2]

    if i <= oy && j <= ox
        out[i, j] = _ka_cubic_sample(a, yval[i], xval[j])
    end
end

@kernel function _ka_rotate_linear_kernel!(
    out,
    old_image,
    c,
    s,
    cx,
    cy,
    sx,
    sy,
    ny::Int,
    nx::Int,
)
    I = @index(Global, NTuple)
    i = I[1]
    j = I[2]

    if i <= ny && j <= nx
        x = j - cx - sx
        y = i - cy - sy
        xr = c * x - s * y + cx
        yr = s * x + c * y + cy
        out[i, j] = _ka_bilinear_sample(old_image, yr, xr)
    end
end

@kernel function _ka_rotate_cubic_kernel!(
    out,
    old_image,
    c,
    s,
    cx,
    cy,
    sx,
    sy,
    ny::Int,
    nx::Int,
)
    I = @index(Global, NTuple)
    i = I[1]
    j = I[2]

    if i <= ny && j <= nx
        x = j - cx - sx
        y = i - cy - sy
        xr = c * x - s * y + cx
        yr = s * x + c * y + cy
        out[i, j] = _ka_cubic_sample(old_image, yr, xr)
    end
end

@inline function ka_apply_shifted_mask!(
    field::AbstractMatrix{<:Complex},
    mask::AbstractMatrix{<:Real};
    invert::Bool=false,
)
    ny, nx = size(field)
    backend = AK.get_backend(field)
    _ka_apply_shifted_mask_kernel!(backend, (16, 16))(field, mask, ny ÷ 2, nx ÷ 2, ny, nx, invert; ndrange=(ny, nx))
    AK.synchronize(backend)
    return field
end

@inline function ka_copy_shifted_complex!(
    out::AbstractMatrix{<:Complex},
    field::AbstractMatrix{<:Complex},
    r0::Int,
    c0::Int,
)
    ny, nx = size(field)
    oy, ox = size(out)
    backend = AK.get_backend(out)
    _ka_copy_shifted_complex_kernel!(backend, (16, 16))(out, field, r0, c0, ny ÷ 2, nx ÷ 2, ny, nx, oy, ox; ndrange=(oy, ox))
    AK.synchronize(backend)
    return out
end

@inline function ka_cubic_conv_grid!(
    out::AbstractMatrix,
    a::AbstractMatrix,
    xval::AbstractVector,
    yval::AbstractVector,
)
    oy, ox = size(out)
    backend = AK.get_backend(out)
    _ka_cubic_conv_grid_kernel!(backend, (16, 16))(out, a, xval, yval, oy, ox; ndrange=(oy, ox))
    AK.synchronize(backend)
    return out
end

@inline function ka_rotate_linear!(
    out::AbstractMatrix,
    old_image::AbstractMatrix,
    c,
    s,
    cx,
    cy,
    sx,
    sy,
)
    ny, nx = size(old_image)
    backend = AK.get_backend(out)
    _ka_rotate_linear_kernel!(backend, (16, 16))(out, old_image, c, s, cx, cy, sx, sy, ny, nx; ndrange=(ny, nx))
    AK.synchronize(backend)
    return out
end

@inline function ka_rotate_cubic!(
    out::AbstractMatrix,
    old_image::AbstractMatrix,
    c,
    s,
    cx,
    cy,
    sx,
    sy,
)
    ny, nx = size(old_image)
    backend = AK.get_backend(out)
    _ka_rotate_cubic_kernel!(backend, (16, 16))(out, old_image, c, s, cx, cy, sx, sy, ny, nx; ndrange=(ny, nx))
    AK.synchronize(backend)
    return out
end

@inline function ka_copy_shifted_intensity!(
    out::AbstractMatrix,
    field::AbstractMatrix{<:Complex},
    r0::Int,
    c0::Int,
)
    ny, nx = size(field)
    oy, ox = size(out)
    backend = AK.get_backend(out)
    _ka_copy_shifted_intensity_kernel!(backend, (16, 16))(out, field, r0, c0, ny ÷ 2, nx ÷ 2, ny, nx, oy, ox; ndrange=(oy, ox))
    AK.synchronize(backend)
    return out
end
