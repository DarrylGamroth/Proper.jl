struct RectangleOptions{T<:AbstractFloat}
    norm::Bool
    dark::Bool
    rotation::T
end

@inline function RectangleOptions(kwargs::Base.Iterators.Pairs)
    return RectangleOptions{Float64}(
        kw_lookup_bool(kwargs, :NORM, false),
        kw_lookup_bool(kwargs, :DARK, false),
        kw_lookup_float(kwargs, :ROTATION, 0.0),
    )
end

function _prop_rectangle!(
    ::GeometryLoopExecStyle,
    image::AbstractMatrix,
    wf::WaveFront,
    xsize::Real,
    ysize::Real,
    xc::Real,
    yc::Real,
    opts::RectangleOptions,
)
    n = prop_get_gridsize(wf)
    size(image) == (n, n) || throw(ArgumentError("output size must match wavefront grid"))
    T = eltype(image)
    dx = T(prop_get_sampling(wf))
    beamrad = T(prop_get_beamradius(wf))
    pr = beamrad / dx

    θ = T(deg2rad(opts.rotation))
    sθ, cθ = sincos(θ)

    xcp = T(n ÷ 2) + (opts.norm ? T(xc) * pr : T(xc) / dx)
    ycp = T(n ÷ 2) + (opts.norm ? T(yc) * pr : T(yc) / dx)
    xrp = (opts.norm ? T(xsize) * pr : T(xsize) / dx) / T(2)
    yrp = (opts.norm ? T(ysize) * pr : T(ysize) / dx) / T(2)

    fill!(image, zero(T))

    # Bounding box from rotated corners.
    xp = (-xrp, -xrp, xrp, xrp)
    yp = (-yrp, yrp, yrp, -yrp)
    xbox = ntuple(i -> xp[i] * cθ - yp[i] * sθ + xcp, 4)
    ybox = ntuple(i -> xp[i] * sθ + yp[i] * cθ + ycp, 4)

    minx = max(0, floor(Int, minimum(xbox) - 1))
    maxx = min(n - 1, ceil(Int, maximum(xbox) + 1))
    miny = max(0, floor(Int, minimum(ybox) - 1))
    maxy = min(n - 1, ceil(Int, maximum(ybox) + 1))

    nsub = antialias_subsampling()
    inv_sub = inv(T(nsub * nsub))

    @inbounds for ypix in miny:maxy
        y0 = ypix - ycp
        for xpix in minx:maxx
            x0 = xpix - xcp
            cnt = 0
            for ys in 1:nsub
                yo = y0 + (ys - (nsub + 1) / 2) / nsub
                for xs in 1:nsub
                    xo = x0 + (xs - (nsub + 1) / 2) / nsub
                    xr = xo * cθ - yo * sθ
                    yr = xo * sθ + yo * cθ
                    cnt += (abs(xr) <= xrp && abs(yr) <= yrp)
                end
            end
            image[ypix + 1, xpix + 1] = cnt * inv_sub
        end
    end

    if opts.dark
        image .= 1 .- image
    end

    return image
end

function _prop_rectangle!(
    ::GeometryKAExecStyle,
    image::AbstractMatrix,
    wf::WaveFront,
    xsize::Real,
    ysize::Real,
    xc::Real,
    yc::Real,
    opts::RectangleOptions,
)
    n = prop_get_gridsize(wf)
    size(image) == (n, n) || throw(ArgumentError("output size must match wavefront grid"))
    T = eltype(image)
    dx = T(prop_get_sampling(wf))
    beamrad = T(prop_get_beamradius(wf))
    pr = beamrad / dx

    θ = T(deg2rad(opts.rotation))
    sθ, cθ = sincos(θ)

    xcp = T(n ÷ 2) + (opts.norm ? T(xc) * pr : T(xc) / dx)
    ycp = T(n ÷ 2) + (opts.norm ? T(yc) * pr : T(yc) / dx)
    xrp = (opts.norm ? T(xsize) * pr : T(xsize) / dx) / T(2)
    yrp = (opts.norm ? T(ysize) * pr : T(ysize) / dx) / T(2)

    xp = (-xrp, -xrp, xrp, xrp)
    yp = (-yrp, yrp, yrp, -yrp)
    xbox = ntuple(i -> xp[i] * cθ - yp[i] * sθ + xcp, 4)
    ybox = ntuple(i -> xp[i] * sθ + yp[i] * cθ + ycp, 4)

    minx = max(0, floor(Int, minimum(xbox) - one(T)))
    maxx = min(n - 1, ceil(Int, maximum(xbox) + one(T)))
    miny = max(0, floor(Int, minimum(ybox) - one(T)))
    maxy = min(n - 1, ceil(Int, maximum(ybox) + one(T)))

    return ka_rectangle_mask!(
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
        maxy;
        dark=opts.dark,
        nsub=antialias_subsampling(),
    )
end

function _prop_rectangle!(
    image::AbstractMatrix,
    wf::WaveFront,
    xsize::Real,
    ysize::Real,
    xc::Real,
    yc::Real,
    opts::RectangleOptions,
)
    return _prop_rectangle!(geometry_exec_style(typeof(image), size(image, 1), size(image, 2)), image, wf, xsize, ysize, xc, yc, opts)
end

function _prop_rectangle(
    wf::WaveFront,
    xsize::Real,
    ysize::Real,
    xc::Real,
    yc::Real,
    opts::RectangleOptions,
)
    RT = real(eltype(wf.field))
    image = similar(wf.field, RT, prop_get_gridsize(wf), prop_get_gridsize(wf))
    return _prop_rectangle!(image, wf, xsize, ysize, xc, yc, opts)
end

"""
    prop_rectangle!(image, wf, xsize, ysize, xc=0, yc=0; kwargs...)

Write an antialiased filled rectangle mask into `image`.

# Arguments
- `image`: destination mask array
- `wf`: wavefront defining the grid and sampling
- `xsize`, `ysize`: rectangle size
- `xc`, `yc`: rectangle center

# Keywords
- `NORM`: interpret size and center in normalized beam-radius units
- `DARK`: return the complementary dark mask
- `ROTATION`: rectangle rotation angle in degrees
"""
function prop_rectangle!(
    image::AbstractMatrix,
    wf::WaveFront,
    xsize::Real,
    ysize::Real,
    xc::Real=0.0,
    yc::Real=0.0;
    kwargs...,
)
    return _prop_rectangle!(image, wf, xsize, ysize, xc, yc, RectangleOptions(kwargs))
end

"""Return an antialiased filled rectangle mask."""
function prop_rectangle(
    wf::WaveFront,
    xsize::Real,
    ysize::Real,
    xc::Real=0.0,
    yc::Real=0.0;
    kwargs...,
)
    opts = RectangleOptions(kwargs)
    return _prop_rectangle(wf, xsize, ysize, xc, yc, opts)
end
