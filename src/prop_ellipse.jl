struct EllipseOptions{T<:AbstractFloat}
    norm::Bool
    dark::Bool
    rotation::T
end

struct EllipseGeometry{T<:AbstractFloat}
    xcenter_pix::T
    ycenter_pix::T
    xrad_pix::T
    yrad_pix::T
    sint::T
    cost::T
    threshold_hi::T
    threshold_lo::T
    limit::T
    minx_pix::Int
    maxx_pix::Int
    miny_pix::Int
    maxy_pix::Int
end

@inline function EllipseOptions(::Type{T}, kwargs::Base.Iterators.Pairs) where {T<:AbstractFloat}
    return EllipseOptions{T}(
        kw_lookup_bool(kwargs, :NORM, false),
        kw_lookup_bool(kwargs, :DARK, false),
        T(kw_lookup_float(kwargs, :ROTATION, 0.0)),
    )
end

@inline EllipseOptions(kwargs::Base.Iterators.Pairs) = EllipseOptions(Float64, kwargs)

@inline function ellipse_geometry(
    ::Type{T},
    wf::WaveFront,
    xradius::Real,
    yradius::Real,
    xc::Real,
    yc::Real,
    opts::EllipseOptions,
) where {T<:AbstractFloat}
    n = prop_get_gridsize(wf)
    dx = T(prop_get_sampling(wf))
    beamrad_pix = T(prop_get_beamradius(wf)) / dx

    xcenter_pix = T(n ÷ 2)
    ycenter_pix = T(n ÷ 2)
    xrad_pix = zero(T)
    yrad_pix = zero(T)
    if opts.norm
        xcenter_pix += T(xc) * beamrad_pix
        ycenter_pix += T(yc) * beamrad_pix
        xrad_pix = T(xradius) * beamrad_pix
        yrad_pix = T(yradius) * beamrad_pix
    else
        xcenter_pix += T(xc) / dx
        ycenter_pix += T(yc) / dx
        xrad_pix = T(xradius) / dx
        yrad_pix = T(yradius) / dx
    end

    t = T(deg2rad(opts.rotation))
    sint = sin(t)
    cost = cos(t)

    minx = typemax(T)
    maxx = typemin(T)
    miny = typemax(T)
    maxy = typemin(T)
    @inbounds for (xp, yp) in ((-xrad_pix, -yrad_pix), (xrad_pix, -yrad_pix), (xrad_pix, yrad_pix), (-xrad_pix, yrad_pix))
        xr = xp * cost - yp * sint + xcenter_pix
        yr = xp * sint + yp * cost + ycenter_pix
        minx = min(minx, xr)
        maxx = max(maxx, xr)
        miny = min(miny, yr)
        maxy = max(maxy, yr)
    end

    minx_pix = clamp(round(Int, minx) - 1, 0, n - 1)
    maxx_pix = clamp(round(Int, maxx) + 1, 0, n - 1)
    miny_pix = clamp(round(Int, miny) - 1, 0, n - 1)
    maxy_pix = clamp(round(Int, maxy) + 1, 0, n - 1)

    delx = inv(xrad_pix)
    dely = inv(yrad_pix)
    drx = delx * cost - dely * sint
    dry = delx * sint + dely * cost
    dr = max(abs(drx), abs(dry))

    return EllipseGeometry(
        xcenter_pix,
        ycenter_pix,
        xrad_pix,
        yrad_pix,
        sint,
        cost,
        one(T) + dr,
        one(T) - dr,
        T(1 + 1e-10),
        minx_pix,
        maxx_pix,
        miny_pix,
        maxy_pix,
    )
end

@inline function ellipse_mask_fraction(x0::T, y0::T, geom::EllipseGeometry{T}, nsub::Int) where {T<:AbstractFloat}
    xr = (x0 * geom.cost - y0 * geom.sint) / geom.xrad_pix
    yr = (x0 * geom.sint + y0 * geom.cost) / geom.yrad_pix
    rv = sqrt(xr * xr + yr * yr)

    if rv > geom.threshold_hi
        return zero(T)
    elseif rv <= geom.threshold_lo
        return one(T)
    end

    cnt = 0
    @inbounds for oy_i in 1:nsub
        ys = y0 + _ka_subsample_offset(oy_i, nsub, one(T))
        for ox_i in 1:nsub
            xs = x0 + _ka_subsample_offset(ox_i, nsub, one(T))
            xsv = (xs * geom.cost - ys * geom.sint) / geom.xrad_pix
            ysv = (xs * geom.sint + ys * geom.cost) / geom.yrad_pix
            cnt += ((xsv * xsv + ysv * ysv) <= geom.limit)
        end
    end

    return T(cnt) / T(nsub * nsub)
end

@inline function ellipse_mask_value(x0::T, y0::T, geom::EllipseGeometry{T}, dark::Bool, nsub::Int) where {T<:AbstractFloat}
    pixval = ellipse_mask_fraction(x0, y0, geom, nsub)
    return dark ? (one(T) - pixval) : pixval
end

abstract type ShiftedEllipseApplyExecStyle end
struct ShiftedEllipseMaskExecStyle <: ShiftedEllipseApplyExecStyle end
struct ShiftedEllipseKAExecStyle <: ShiftedEllipseApplyExecStyle end

@inline shifted_ellipse_apply_exec_style(::GeometryKAStyle) = ShiftedEllipseKAExecStyle()
@inline shifted_ellipse_apply_exec_style(::GeometryKernelStyle) = ShiftedEllipseMaskExecStyle()
@inline shifted_ellipse_apply_exec_style(::Type{A}) where {A<:AbstractArray} = shifted_ellipse_apply_exec_style(geometry_kernel_style(A))

function _prop_ellipse!(
    ::GeometryLoopExecStyle,
    image::AbstractMatrix,
    wf::WaveFront,
    xradius::Real,
    yradius::Real,
    xc::Real,
    yc::Real,
    opts::EllipseOptions,
)
    n = prop_get_gridsize(wf)
    size(image) == (n, n) || throw(ArgumentError("output size must match wavefront grid"))

    nsub = antialias_subsampling()
    T = eltype(image)
    geom = ellipse_geometry(T, wf, xradius, yradius, xc, yc, opts)
    fill!(image, opts.dark ? one(T) : zero(T))

    @inbounds for ypix in geom.miny_pix:geom.maxy_pix
        y0 = T(ypix) - geom.ycenter_pix
        for xpix in geom.minx_pix:geom.maxx_pix
            x0 = T(xpix) - geom.xcenter_pix
            image[ypix + 1, xpix + 1] = ellipse_mask_value(x0, y0, geom, opts.dark, nsub)
        end
    end

    return image
end

function _prop_ellipse!(
    ::GeometryKAExecStyle,
    image::AbstractMatrix,
    wf::WaveFront,
    xradius::Real,
    yradius::Real,
    xc::Real,
    yc::Real,
    opts::EllipseOptions,
)
    n = prop_get_gridsize(wf)
    size(image) == (n, n) || throw(ArgumentError("output size must match wavefront grid"))

    T = eltype(image)
    geom = ellipse_geometry(T, wf, xradius, yradius, xc, yc, opts)

    return ka_ellipse_mask!(
        image,
        geom.xcenter_pix,
        geom.ycenter_pix,
        geom.xrad_pix,
        geom.yrad_pix,
        geom.sint,
        geom.cost,
        geom.threshold_hi,
        geom.threshold_lo,
        geom.limit;
        dark=opts.dark,
        nsub=antialias_subsampling(),
    )
end

function _prop_ellipse!(
    image::AbstractMatrix,
    wf::WaveFront,
    xradius::Real,
    yradius::Real,
    xc::Real,
    yc::Real,
    opts::EllipseOptions,
)
    return _prop_ellipse!(geometry_exec_style(typeof(image), size(image, 1), size(image, 2)), image, wf, xradius, yradius, xc, yc, opts)
end

function _prop_ellipse(
    wf::WaveFront,
    xradius::Real,
    yradius::Real,
    xc::Real,
    yc::Real,
    opts::EllipseOptions,
)
    RT = real(eltype(wf.field))
    image = similar(wf.field, RT, prop_get_gridsize(wf), prop_get_gridsize(wf))
    return _prop_ellipse!(image, wf, xradius, yradius, xc, yc, opts)
end

function prop_ellipse!(
    image::AbstractMatrix,
    wf::WaveFront,
    xradius::Real,
    yradius::Real,
    xc::Real=0.0,
    yc::Real=0.0;
    kwargs...,
)
    opts = EllipseOptions(eltype(image), kwargs)
    return _prop_ellipse!(image, wf, xradius, yradius, xc, yc, opts)
end

"""Creates an image containing an antialiased filled ellipse."""
function prop_ellipse(
    wf::WaveFront,
    xradius::Real,
    yradius::Real,
    xc::Real=0.0,
    yc::Real=0.0;
    kwargs...,
)
    opts = EllipseOptions(real(eltype(wf.field)), kwargs)
    return _prop_ellipse(wf, xradius, yradius, xc, yc, opts)
end

@inline function _apply_shifted_ellipse_loop!(
    field::AbstractMatrix{Complex{T}},
    geom::EllipseGeometry{T},
    dark::Bool,
    invert::Bool,
    nsub::Int,
) where {T<:AbstractFloat}
    ny, nx = size(field)
    sy = ny ÷ 2
    sx = nx ÷ 2

    @inbounds for j in 1:nx
        js = mod1(j + sx, nx)
        x0 = T(js - 1) - geom.xcenter_pix
        for i in 1:ny
            is = mod1(i + sy, ny)
            y0 = T(is - 1) - geom.ycenter_pix
            maskval = ellipse_mask_value(x0, y0, geom, dark, nsub)
            field[i, j] *= invert ? (one(T) - maskval) : maskval
        end
    end

    return field
end

@inline function _apply_shifted_ellipse!(
    ::ShiftedEllipseMaskExecStyle,
    wf::WaveFront,
    xradius::Real,
    yradius::Real,
    xc::Real,
    yc::Real,
    opts::EllipseOptions,
    invert::Bool,
)
    ny, nx = size(wf.field)
    m = ensure_mask_buffer!(wf.workspace.mask, ny, nx)
    _prop_ellipse!(m, wf, xradius, yradius, xc, yc, opts)
    _apply_shifted_mask!(wf.field, m; invert=invert)
    return wf
end

@inline function _apply_shifted_ellipse!(
    ::ShiftedEllipseKAExecStyle,
    wf::WaveFront{T},
    xradius::Real,
    yradius::Real,
    xc::Real,
    yc::Real,
    opts::EllipseOptions,
    invert::Bool,
) where {T<:AbstractFloat}
    geom = ellipse_geometry(T, wf, xradius, yradius, xc, yc, opts)
    ka_apply_shifted_ellipse!(
        wf.field,
        geom.xcenter_pix,
        geom.ycenter_pix,
        geom.xrad_pix,
        geom.yrad_pix,
        geom.sint,
        geom.cost,
        geom.threshold_hi,
        geom.threshold_lo,
        geom.limit;
        dark=opts.dark,
        invert=invert,
        nsub=antialias_subsampling(),
    )
    return wf
end

@inline function _apply_shifted_ellipse!(
    wf::WaveFront,
    xradius::Real,
    yradius::Real,
    xc::Real,
    yc::Real,
    opts::EllipseOptions,
    invert::Bool,
)
    sty = shifted_ellipse_apply_exec_style(typeof(wf.field))
    return _apply_shifted_ellipse!(sty, wf, xradius, yradius, xc, yc, opts, invert)
end
