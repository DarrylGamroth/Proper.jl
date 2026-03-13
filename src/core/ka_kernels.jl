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

@inline function _ka_shifted_index_0based(p::Int, n::Int)
    c = n ÷ 2
    cut = n - c - 1
    return ifelse(p <= cut, p, p - n)
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

@inline function _ka_szoom_round(x::T) where {T<:AbstractFloat}
    return x < zero(T) ? floor(x) : ceil(x)
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
    @Const(mask),
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

@kernel function _ka_apply_shifted_ellipse_kernel!(
    field,
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
    invert::Bool,
    minx_pix::Int,
    maxx_pix::Int,
    miny_pix::Int,
    maxy_pix::Int,
    sy::Int,
    sx::Int,
    ny::Int,
    nx::Int,
)
    I = @index(Global, NTuple)
    i = I[1]
    j = I[2]

    if i <= ny && j <= nx
        @uniform begin
            T = typeof(xcenter_pix)
            inv_nsub2 = inv(T(nsub * nsub))
            outside_factor = xor(dark, invert) ? one(T) : zero(T)
            inside_factor = xor(dark, invert) ? zero(T) : one(T)
        end
        is = i + sy
        js = j + sx
        is = ifelse(is > ny, is - ny, is)
        js = ifelse(js > nx, js - nx, js)
        is0 = is - 1
        js0 = js - 1

        factor = outside_factor

        if !(js0 < minx_pix || js0 > maxx_pix || is0 < miny_pix || is0 > maxy_pix)
            x0 = T(js - 1) - xcenter_pix
            y0 = T(is - 1) - ycenter_pix

            xr = (x0 * cost - y0 * sint) / xrad_pix
            yr = (x0 * sint + y0 * cost) / yrad_pix
            rv = sqrt(xr * xr + yr * yr)

            if rv <= threshold_lo
                factor = inside_factor
            elseif rv <= threshold_hi
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
                pixval = T(cnt) * inv_nsub2
                maskval = dark ? (one(T) - pixval) : pixval
                factor = invert ? (one(T) - maskval) : maskval
            end
        end
        if factor == zero(T)
            field[i, j] = zero(eltype(field))
        elseif factor != one(T)
            field[i, j] *= factor
        end
    end
end

@kernel function _ka_apply_shifted_circle_kernel!(
    field,
    xoffset_pix,
    yoffset_pix,
    rad_pix,
    threshold_hi2,
    threshold_lo2,
    limit2,
    nsub::Int,
    dark::Bool,
    invert::Bool,
    ny::Int,
    nx::Int,
)
    I = @index(Global, NTuple)
    i = I[1]
    j = I[2]

    if i <= ny && j <= nx
        @uniform begin
            T = typeof(rad_pix)
            inv_nsub2 = inv(T(nsub * nsub))
            outside_factor = xor(dark, invert) ? one(T) : zero(T)
            inside_factor = xor(dark, invert) ? zero(T) : one(T)
        end
        x0 = T(_ka_shifted_index_0based(j - 1, nx)) - xoffset_pix
        y0 = T(_ka_shifted_index_0based(i - 1, ny)) - yoffset_pix
        r2 = x0 * x0 + y0 * y0

        factor = outside_factor
        if r2 <= threshold_lo2
            factor = inside_factor
        elseif r2 <= threshold_hi2
            cnt = 0
            for oy_i in 1:nsub
                ys = y0 + _ka_subsample_offset(oy_i, nsub, one(T))
                for ox_i in 1:nsub
                    xs = x0 + _ka_subsample_offset(ox_i, nsub, one(T))
                    cnt += ((xs * xs + ys * ys) <= limit2)
                end
            end
            pixval = T(cnt) * inv_nsub2
            maskval = dark ? (one(T) - pixval) : pixval
            factor = invert ? (one(T) - maskval) : maskval
        end

        if factor == zero(T)
            field[i, j] = zero(eltype(field))
        elseif factor != one(T)
            field[i, j] *= factor
        end
    end
end

@kernel function _ka_copy_shifted_complex_kernel!(
    out,
    @Const(field),
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
    @Const(field),
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

@kernel function _ka_apply_qphase_kernel!(
    field,
    k,
    dx,
    ny::Int,
    nx::Int,
)
    I = @index(Global, NTuple)
    i = I[1]
    j = I[2]

    if i <= ny && j <= nx
        T = typeof(dx)
        x = T(_ka_shifted_index_0based(j - 1, nx)) * dx
        y = T(_ka_shifted_index_0based(i - 1, ny)) * dx
        field[i, j] *= cis(k * (x * x + y * y))
    end
end

@kernel function _ka_apply_frequency_phase_kernel!(
    field,
    kphase,
    inv_dx_y,
    inv_dx_x,
    ny::Int,
    nx::Int,
)
    I = @index(Global, NTuple)
    i = I[1]
    j = I[2]

    if i <= ny && j <= nx
        T = typeof(inv_dx_y)
        fy = T(_ka_shifted_index_0based(i - 1, ny)) * inv_dx_y
        fx = T(_ka_shifted_index_0based(j - 1, nx)) * inv_dx_x
        field[i, j] *= cis(kphase * (fx * fx + fy * fy))
    end
end

@kernel function _ka_fill_affine_axis_kernel!(
    out,
    origin,
    scale,
    offset,
    n::Int,
)
    i = @index(Global, Linear)
    if i <= n
        T = typeof(scale)
        out[i] = (T(i - 1) - origin) * scale + offset
    end
end

@kernel function _ka_cubic_conv_grid_kernel!(
    out,
    @Const(a),
    @Const(xval),
    @Const(yval),
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
    @Const(old_image),
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
    @Const(old_image),
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
    ystart::Int,
    xstart::Int,
    boxny::Int,
    boxnx::Int,
    nsub::Int,
)
    I = @index(Global, NTuple)
    bi = I[1]
    bj = I[2]

    if bi <= boxny && bj <= boxnx
        @uniform begin
            T = typeof(xcp)
            inv_nsub2 = inv(T(nsub * nsub))
        end
        i = ystart + bi - 1
        j = xstart + bj - 1
        xpix = j - 1
        ypix = i - 1
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

        image[i, j] = T(cnt) * inv_nsub2
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
        @uniform begin
            T = typeof(xcenter_pix)
            inv_nsub2 = inv(T(nsub * nsub))
        end
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
            T(cnt) * inv_nsub2
        end

        image[i, j] = dark ? (one(T) - pixval) : pixval
    end
end

@kernel function _ka_irregular_polygon_mask_kernel!(
    image,
    @Const(xv),
    @Const(yv),
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
        @uniform begin
            T = typeof(dx)
            inv_nsub2 = inv(T(nsub * nsub))
        end
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

        pixval = T(cnt) * inv_nsub2
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
    ystart::Int,
    xstart::Int,
    boxny::Int,
    boxnx::Int,
)
    I = @index(Global, NTuple)
    bi = I[1]
    bj = I[2]

    if bi <= boxny && bj <= boxnx
        T = typeof(dx)
        i = ystart + bi - 1
        j = xstart + bj - 1
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

@kernel function _ka_szoom_table_kernel!(
    table,
    mag,
    n_out::Int,
    k::Int,
    dk::Int,
)
    I = @index(Global, NTuple)
    row = I[1]
    idx = I[2]

    if row <= n_out && idx <= k
        T = typeof(mag)
        xin = T(row - 1 - (n_out ÷ 2)) / mag
        xphase = xin - _ka_szoom_round(xin)
        x = T(idx - 1 - (k ÷ 2)) - xphase

        if abs(x) <= T(dk)
            xpi = x * T(pi)
            if xpi == zero(T)
                table[row, idx] = one(T)
            else
                table[row, idx] = (sin(xpi) / xpi) * (sin(xpi / T(dk)) / (xpi / T(dk)))
            end
        else
            table[row, idx] = zero(T)
        end
    end
end

@kernel function _ka_szoom_apply_kernel!(
    out,
    @Const(image_in),
    @Const(table),
    mag,
    n_in::Int,
    n_out::Int,
    k::Int,
)
    I = @index(Global, NTuple)
    row_out = I[1]
    col_out = I[2]

    if row_out <= n_out && col_out <= n_out
        T = typeof(float(real(zero(eltype(out)))))
        yin = T(row_out - 1 - (n_out ÷ 2)) / mag
        ypix = _ka_szoom_round(yin) + T(n_in ÷ 2)
        y1 = Int(ypix) - (k ÷ 2)
        y2_excl = Int(ypix) + (k ÷ 2) + 1

        xin = T(col_out - 1 - (n_out ÷ 2)) / mag
        xpix = _ka_szoom_round(xin) + T(n_in ÷ 2)
        x1 = Int(xpix) - (k ÷ 2)
        x2_excl = Int(xpix) + (k ÷ 2) + 1

        if y1 < 0 || y2_excl > n_in || x1 < 0 || x2_excl > n_in
            out[row_out, col_out] = zero(eltype(out))
        elseif eltype(image_in) <: Complex
            acc_re = zero(T)
            acc_im = zero(T)
            for co in 0:(k - 1)
                col = x1 + co + 1
                s_re = zero(T)
                s_im = zero(T)
                for ro in 0:(k - 1)
                    row = y1 + ro + 1
                    wrow = table[row_out, ro + 1]
                    z = image_in[row, col]
                    s_re += T(real(z)) * wrow
                    s_im += T(imag(z)) * wrow
                end
                wcol = table[col_out, co + 1]
                acc_re += s_re * wcol
                acc_im += s_im * wcol
            end
            out[row_out, col_out] = complex(acc_re, acc_im)
        else
            acc = zero(T)
            for co in 0:(k - 1)
                col = x1 + co + 1
                scol = zero(T)
                for ro in 0:(k - 1)
                    row = y1 + ro + 1
                    scol += T(image_in[row, col]) * table[row_out, ro + 1]
                end
                acc += scol * table[col_out, co + 1]
            end
            out[row_out, col_out] = acc
        end
    end
end

@kernel function _ka_pixellate_kernel!(
    out,
    @Const(img),
    f::Int,
    scale,
    ny2::Int,
    nx2::Int,
)
    I = @index(Global, NTuple)
    i = I[1]
    j = I[2]

    if i <= ny2 && j <= nx2
        ys = (i - 1) * f + 1
        xs = (j - 1) * f + 1
        acc = zero(eltype(out))
        for xoff in 0:(f - 1)
            for yoff in 0:(f - 1)
                acc += img[ys + yoff, xs + xoff]
            end
        end
        out[i, j] = acc * scale
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
    return field
end

@inline function ka_apply_shifted_ellipse!(
    field::AbstractMatrix{<:Complex},
    xcenter_pix,
    ycenter_pix,
    xrad_pix,
    yrad_pix,
    sint,
    cost,
    threshold_hi,
    threshold_lo,
    limit;
    minx_pix::Int=0,
    maxx_pix::Int=size(field, 2) - 1,
    miny_pix::Int=0,
    maxy_pix::Int=size(field, 1) - 1,
    dark::Bool=false,
    invert::Bool=false,
    nsub::Int=1,
)
    ny, nx = size(field)
    backend = AK.get_backend(field)
    _ka_apply_shifted_ellipse_kernel!(
        backend,
        (16, 16),
    )(
        field,
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
        invert,
        minx_pix,
        maxx_pix,
        miny_pix,
        maxy_pix,
        ny ÷ 2,
        nx ÷ 2,
        ny,
        nx;
        ndrange=(ny, nx),
    )
    return field
end

@inline function ka_apply_shifted_circle!(
    field::AbstractMatrix{<:Complex},
    xoffset_pix,
    yoffset_pix,
    rad_pix,
    threshold_hi2,
    threshold_lo2,
    limit2;
    dark::Bool=false,
    invert::Bool=false,
    nsub::Int=1,
)
    ny, nx = size(field)
    backend = AK.get_backend(field)
    _ka_apply_shifted_circle_kernel!(
        backend,
        (16, 16),
    )(
        field,
        xoffset_pix,
        yoffset_pix,
        rad_pix,
        threshold_hi2,
        threshold_lo2,
        limit2,
        nsub,
        dark,
        invert,
        ny,
        nx;
        ndrange=(ny, nx),
    )
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
    return out
end

@inline function ka_apply_qphase!(
    field::AbstractMatrix{<:Complex},
    k,
    dx,
)
    ny, nx = size(field)
    backend = AK.get_backend(field)
    _ka_apply_qphase_kernel!(backend, (16, 16))(field, k, dx, ny, nx; ndrange=(ny, nx))
    return field
end

@inline function ka_apply_frequency_phase!(
    field::AbstractMatrix{<:Complex},
    kphase,
    dx,
)
    ny, nx = size(field)
    T = typeof(dx)
    inv_dx_y = inv(T(ny) * T(dx))
    inv_dx_x = inv(T(nx) * T(dx))
    backend = AK.get_backend(field)
    _ka_apply_frequency_phase_kernel!(backend, (16, 16))(field, kphase, inv_dx_y, inv_dx_x, ny, nx; ndrange=(ny, nx))
    return field
end

@inline function ka_fill_affine_axis!(
    out::AbstractVector,
    origin,
    scale,
    offset,
)
    n = length(out)
    backend = AK.get_backend(out)
    _ka_fill_affine_axis_kernel!(backend, 256)(out, origin, scale, offset, n; ndrange=n)
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
    fill!(image, dark ? one(T) : zero(T))
    boxny = max(maxy - miny + 1, 0)
    boxnx = max(maxx - minx + 1, 0)
    boxny == 0 && return image
    boxnx == 0 && return image
    backend = AK.get_backend(image)
    _ka_rectangle_mask_kernel!(backend, (16, 16))(
        image,
        xcp,
        ycp,
        xrp,
        yrp,
        cθ,
        sθ,
        miny + 1,
        minx + 1,
        boxny,
        boxnx,
        nsub,
        ndrange=(boxny, boxnx),
    )
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
    cx_pix = T(nx ÷ 2) + xc / dx
    cy_pix = T(ny ÷ 2) + yc / dx
    halfw_pix = hw / dx
    halfh_pix = hh / dx
    minx = max(0, floor(Int, cx_pix - halfw_pix - one(T)))
    maxx = min(nx - 1, ceil(Int, cx_pix + halfw_pix + one(T)))
    miny = max(0, floor(Int, cy_pix - halfh_pix - one(T)))
    maxy = min(ny - 1, ceil(Int, cy_pix + halfh_pix + one(T)))
    boxny = max(maxy - miny + 1, 0)
    boxnx = max(maxx - minx + 1, 0)
    fill!(image, zero(T))
    boxny == 0 && return image
    boxnx == 0 && return image
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
        miny + 1,
        minx + 1,
        boxny,
        boxnx;
        ndrange=(boxny, boxnx),
    )
    return image
end

@inline function ka_szoom_table!(
    table::AbstractMatrix{T},
    mag::T,
    n_out::Int,
    k::Int,
    dk::Int,
) where {T<:AbstractFloat}
    backend = AK.get_backend(table)
    _ka_szoom_table_kernel!(backend, (16, 16))(table, mag, n_out, k, dk; ndrange=size(table))
    return table
end

@inline function ka_szoom_apply!(
    out::AbstractMatrix,
    image_in::AbstractMatrix,
    table::AbstractMatrix,
    mag,
)
    n_out = size(out, 1)
    k = size(table, 2)
    backend = AK.get_backend(out)
    _ka_szoom_apply_kernel!(backend, (16, 16))(out, image_in, table, mag, size(image_in, 1), n_out, k; ndrange=size(out))
    return out
end

@inline function ka_pixellate!(
    out::AbstractMatrix,
    img::AbstractMatrix,
    f::Int,
)
    ny2, nx2 = size(out)
    backend = AK.get_backend(out)
    T = typeof(float(real(zero(eltype(out)))))
    scale = inv(T(f * f))
    _ka_pixellate_kernel!(backend, (16, 16))(out, img, f, scale, ny2, nx2; ndrange=(ny2, nx2))
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
    return out
end
