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

@inline function _ka_subsample_offset(idx::Int, nsub::Int, scale)
    T = typeof(scale)
    return (T(idx) - (T(nsub) + one(T)) / T(2)) / T(nsub) * scale
end

@inline function _ka_point_in_poly(x::Real, y::Real, xv::AbstractVector, yv::AbstractVector)
    inside = false
    n = length(xv)
    j = n
    epsy = eps(promote_type(typeof(x), typeof(y), eltype(xv), eltype(yv)))
    @inbounds for i in 1:n
        yi = yv[i]
        yj = yv[j]
        if (yi > y) != (yj > y)
            xcross = (xv[j] - xv[i]) * (y - yi) / (yj - yi + epsy) + xv[i]
            inside = ifelse(x < xcross, !inside, inside)
        end
        j = i
    end
    return inside
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

@kernel function _ka_rectangle_mask_kernel!(
    image,
    xcp,
    ycp,
    xrp,
    yrp,
    cθ,
    sθ,
    minx::Int,
    maxx::Int,
    miny::Int,
    maxy::Int,
    nsub::Int,
    dark::Bool,
    ny::Int,
    nx::Int,
)
    I = @index(Global, NTuple)
    i = I[1]
    j = I[2]

    if i <= ny && j <= nx
        T = typeof(xcp)
        xpix = j - 1
        ypix = i - 1
        pixval = zero(T)

        if minx <= xpix <= maxx && miny <= ypix <= maxy
            x0 = T(xpix) - xcp
            y0 = T(ypix) - ycp
            cnt = 0

            for ys in 1:nsub
                yo = y0 + _ka_subsample_offset(ys, nsub, one(T))
                for xs in 1:nsub
                    xo = x0 + _ka_subsample_offset(xs, nsub, one(T))
                    xr = xo * cθ - yo * sθ
                    yr = xo * sθ + yo * cθ
                    cnt += (abs(xr) <= xrp && abs(yr) <= yrp)
                end
            end

            pixval = T(cnt) / T(nsub * nsub)
        end

        image[i, j] = dark ? (one(T) - pixval) : pixval
    end
end

@kernel function _ka_ellipse_mask_kernel!(
    image,
    xcenter_pix,
    ycenter_pix,
    xrad_pix,
    yrad_pix,
    sint,
    cost,
    threshold_hi,
    threshold_lo,
    limit,
    nsub::Int,
    dark::Bool,
    ny::Int,
    nx::Int,
)
    I = @index(Global, NTuple)
    i = I[1]
    j = I[2]

    if i <= ny && j <= nx
        T = typeof(xcenter_pix)
        x0 = T(j - 1) - xcenter_pix
        y0 = T(i - 1) - ycenter_pix

        xr = (x0 * cost - y0 * sint) / xrad_pix
        yr = (x0 * sint + y0 * cost) / yrad_pix
        rv = sqrt(xr * xr + yr * yr)

        pixval = if rv > threshold_hi
            zero(T)
        elseif rv <= threshold_lo
            one(T)
        else
            cnt = 0
            for oy_i in 1:nsub
                ys = y0 + _ka_subsample_offset(oy_i, nsub, one(T))
                for ox_i in 1:nsub
                    xs = x0 + _ka_subsample_offset(ox_i, nsub, one(T))
                    xsv = (xs * cost - ys * sint) / xrad_pix
                    ysv = (xs * sint + ys * cost) / yrad_pix
                    cnt += ((xsv * xsv + ysv * ysv) <= limit)
                end
            end
            T(cnt) / T(nsub * nsub)
        end

        image[i, j] = dark ? (one(T) - pixval) : pixval
    end
end

@kernel function _ka_irregular_polygon_mask_kernel!(
    image,
    xv,
    yv,
    cx::Int,
    cy::Int,
    dx,
    nsub::Int,
    dark::Bool,
    ny::Int,
    nx::Int,
)
    I = @index(Global, NTuple)
    i = I[1]
    j = I[2]

    if i <= ny && j <= nx
        T = typeof(dx)
        x0 = T(j - 1 - cx) * dx
        y0 = T(i - 1 - cy) * dx
        cnt = 0

        for oy_i in 1:nsub
            oy = _ka_subsample_offset(oy_i, nsub, dx)
            for ox_i in 1:nsub
                ox = _ka_subsample_offset(ox_i, nsub, dx)
                cnt += _ka_point_in_poly(x0 + ox, y0 + oy, xv, yv)
            end
        end

        pixval = T(cnt) / T(nsub * nsub)
        image[i, j] = dark ? (one(T) - pixval) : pixval
    end
end

@kernel function _ka_rounded_rectangle_mask_kernel!(
    image,
    dx,
    xc,
    yc,
    r,
    hw,
    hh,
    cy::Int,
    cx::Int,
    ny::Int,
    nx::Int,
)
    I = @index(Global, NTuple)
    i = I[1]
    j = I[2]

    if i <= ny && j <= nx
        T = typeof(dx)
        x = T(j - 1 - cx) * dx - xc
        y = T(i - 1 - cy) * dx - yc
        qx = abs(x) - (hw - r)
        qy = abs(y) - (hh - r)
        ax = max(qx, zero(T))
        ay = max(qy, zero(T))
        inside = (qx <= 0 && abs(y) <= hh) || (qy <= 0 && abs(x) <= hw) || (ax * ax + ay * ay <= r * r)
        image[i, j] = inside ? one(T) : zero(T)
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

@inline function ka_rectangle_mask!(
    image::AbstractMatrix{T},
    xcp::T,
    ycp::T,
    xrp::T,
    yrp::T,
    cθ::T,
    sθ::T,
    minx::Int,
    maxx::Int,
    miny::Int,
    maxy::Int;
    dark::Bool=false,
    nsub::Int=1,
) where {T<:AbstractFloat}
    ny, nx = size(image)
    backend = AK.get_backend(image)
    _ka_rectangle_mask_kernel!(backend, (16, 16))(
        image,
        xcp,
        ycp,
        xrp,
        yrp,
        cθ,
        sθ,
        minx,
        maxx,
        miny,
        maxy,
        nsub,
        dark,
        ny,
        nx;
        ndrange=(ny, nx),
    )
    AK.synchronize(backend)
    return image
end

@inline function ka_ellipse_mask!(
    image::AbstractMatrix{T},
    xcenter_pix::T,
    ycenter_pix::T,
    xrad_pix::T,
    yrad_pix::T,
    sint::T,
    cost::T,
    threshold_hi::T,
    threshold_lo::T,
    limit::T;
    dark::Bool=false,
    nsub::Int=1,
) where {T<:AbstractFloat}
    ny, nx = size(image)
    backend = AK.get_backend(image)
    _ka_ellipse_mask_kernel!(backend, (16, 16))(
        image,
        xcenter_pix,
        ycenter_pix,
        xrad_pix,
        yrad_pix,
        sint,
        cost,
        threshold_hi,
        threshold_lo,
        limit,
        nsub,
        dark,
        ny,
        nx;
        ndrange=(ny, nx),
    )
    AK.synchronize(backend)
    return image
end

@inline function ka_irregular_polygon_mask!(
    image::AbstractMatrix{T},
    xv::AbstractVector,
    yv::AbstractVector,
    cx::Int,
    cy::Int,
    dx::T;
    dark::Bool=false,
    nsub::Int=1,
) where {T<:AbstractFloat}
    ny, nx = size(image)
    backend = AK.get_backend(image)
    _ka_irregular_polygon_mask_kernel!(backend, (16, 16))(
        image,
        xv,
        yv,
        cx,
        cy,
        dx,
        nsub,
        dark,
        ny,
        nx;
        ndrange=(ny, nx),
    )
    AK.synchronize(backend)
    return image
end

@inline function ka_rounded_rectangle_mask!(
    image::AbstractMatrix{T},
    dx::T,
    xc::T,
    yc::T,
    r::T,
    hw::T,
    hh::T,
) where {T<:AbstractFloat}
    ny, nx = size(image)
    backend = AK.get_backend(image)
    _ka_rounded_rectangle_mask_kernel!(backend, (16, 16))(
        image,
        dx,
        xc,
        yc,
        r,
        hw,
        hh,
        ny ÷ 2,
        nx ÷ 2,
        ny,
        nx;
        ndrange=(ny, nx),
    )
    AK.synchronize(backend)
    return image
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
