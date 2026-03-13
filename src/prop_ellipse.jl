struct EllipseOptions{T<:AbstractFloat}
    norm::Bool
    dark::Bool
    rotation::T
end

@inline function EllipseOptions(kwargs::Base.Iterators.Pairs)
    return EllipseOptions{Float64}(
        kw_lookup_bool(kwargs, :NORM, false),
        kw_lookup_bool(kwargs, :DARK, false),
        kw_lookup_float(kwargs, :ROTATION, 0.0),
    )
end

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
    dx = T(prop_get_sampling(wf))
    beamrad_pix = T(prop_get_beamradius(wf)) / dx
    norm = opts.norm
    dark = opts.dark
    rotation = opts.rotation

    xf = T
    xcenter_pix = xf(n ÷ 2)
    ycenter_pix = xf(n ÷ 2)

    xrad_pix::T = zero(T)
    yrad_pix::T = zero(T)
    if norm
        xcenter_pix += xf(xc) * xf(beamrad_pix)
        ycenter_pix += xf(yc) * xf(beamrad_pix)
        xrad_pix = xf(xradius) * xf(beamrad_pix)
        yrad_pix = xf(yradius) * xf(beamrad_pix)
    else
        xcenter_pix += xf(xc / dx)
        ycenter_pix += xf(yc / dx)
        xrad_pix = xf(xradius / dx)
        yrad_pix = xf(yradius / dx)
    end

    t = xf(deg2rad(rotation))
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

    fill!(image, dark ? one(T) : zero(T))

    threshold_hi = xf(1) + dr
    threshold_lo = xf(1) - dr
    limit = xf(1 + 1e-10)

    nsubpix = xf(nsub * nsub)
    nsub_inv = inv(xf(nsub))
    nsub_ctr = nsub ÷ 2

    @inbounds for ypix in miny_pix:maxy_pix
        y0 = xf(ypix) - ycenter_pix
        for xpix in minx_pix:maxx_pix
            x0 = xf(xpix) - xcenter_pix

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
                    oy = xf(oy_i - 1 - nsub_ctr) * nsub_inv
                    ys = y0 + oy
                    for ox_i in 1:nsub
                        ox = xf(ox_i - 1 - nsub_ctr) * nsub_inv
                        xs = x0 + ox
                        xsv = (xs * cost - ys * sint) / xrad_pix
                        ysv = (xs * sint + ys * cost) / yrad_pix
                        cnt += ((xsv * xsv + ysv * ysv) <= limit)
                    end
                end
                xf(cnt) / nsubpix
            end

            image[ypix + 1, xpix + 1] = dark ? (one(T) - pixval) : pixval
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

    nsub = antialias_subsampling()
    T = eltype(image)
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

    delx = inv(xrad_pix)
    dely = inv(yrad_pix)
    drx = delx * cost - dely * sint
    dry = delx * sint + dely * cost
    dr = max(abs(drx), abs(dry))

    threshold_hi = one(T) + dr
    threshold_lo = one(T) - dr
    limit = T(1 + 1e-10)

    return ka_ellipse_mask!(
        image,
        xcenter_pix,
        ycenter_pix,
        xrad_pix,
        yrad_pix,
        sint,
        cost,
        threshold_hi,
        threshold_lo,
        limit;
        dark=opts.dark,
        nsub=nsub,
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
    opts = EllipseOptions(kwargs)
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
    opts = EllipseOptions(kwargs)
    return _prop_ellipse(wf, xradius, yradius, xc, yc, opts)
end
