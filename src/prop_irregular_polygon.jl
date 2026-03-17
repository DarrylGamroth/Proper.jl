@inline function _point_in_poly(x::Real, y::Real, xv::AbstractVector, yv::AbstractVector)
    inside = false
    n = length(xv)
    j = n
    epsy = eps(promote_type(typeof(y), eltype(yv)))
    @inbounds for i in 1:n
        yi = yv[i]
        yj = yv[j]
        if ((yi > y) != (yj > y))
            xcross = (xv[j] - xv[i]) * (y - yi) / (yj - yi + epsy) + xv[i]
            inside = (x < xcross) ? !inside : inside
        end
        j = i
    end
    return inside
end

struct IrregularPolygonOptions
    dark::Bool
    norm::Bool
end

@inline function IrregularPolygonOptions(kwargs::Base.Iterators.Pairs)
    return IrregularPolygonOptions(
        kw_lookup_bool(kwargs, :DARK, false),
        kw_lookup_bool(kwargs, :NORM, false),
    )
end

function _prop_irregular_polygon(
    ::GeometryLoopExecStyle,
    image::AbstractMatrix,
    wf::WaveFront,
    xverts::AbstractVector,
    yverts::AbstractVector,
    opts::IrregularPolygonOptions,
)
    length(xverts) == length(yverts) || throw(ArgumentError("vertex arrays must have same length"))
    nverts = length(xverts)

    RT = eltype(image)
    beamr = prop_get_beamradius(wf)
    xv, yv = ensure_mask_vertices!(wf.workspace.mask, nverts)
    @inbounds for k in 1:nverts
        xv0 = RT(xverts[k])
        yv0 = RT(yverts[k])
        if opts.norm
            xv[k] = xv0 * RT(beamr)
            yv[k] = yv0 * RT(beamr)
        else
            xv[k] = xv0
            yv[k] = yv0
        end
    end

    ny, nx = size(wf.field)
    size(image) == (ny, nx) || throw(ArgumentError("output size must match wavefront"))
    dx = RT(wf.sampling_m)
    cx = nx ÷ 2
    cy = ny ÷ 2

    subs = antialias_subsampling()
    inv_subs = inv(RT(subs))
    halfsub = RT(subs + 1) / RT(2)
    nsubpix = RT(subs * subs)

    fill!(image, zero(eltype(image)))
    @inbounds for j in 1:nx
        x0 = RT(j - 1 - cx) * dx
        for i in 1:ny
            y0 = RT(i - 1 - cy) * dx
            cnt = 0
            for oy_i in 1:subs
                oy = (RT(oy_i) - halfsub) * inv_subs * dx
                for ox_i in 1:subs
                    ox = (RT(ox_i) - halfsub) * inv_subs * dx
                cnt += _point_in_poly(x0 + ox, y0 + oy, xv, yv) ? 1 : 0
                end
            end
            image[i, j] = RT(cnt) / nsubpix
        end
    end

    if opts.dark
        image .= 1 .- image
    end
    return image
end

function _prop_irregular_polygon(
    ::GeometryKAExecStyle,
    image::AbstractMatrix,
    wf::WaveFront,
    xverts::AbstractVector,
    yverts::AbstractVector,
    opts::IrregularPolygonOptions,
)
    length(xverts) == length(yverts) || throw(ArgumentError("vertex arrays must have same length"))
    nverts = length(xverts)

    RT = eltype(image)
    beamr = prop_get_beamradius(wf)
    xv, yv = ensure_mask_vertices!(wf.workspace.mask, nverts)
    @inbounds for k in 1:nverts
        xv0 = RT(xverts[k])
        yv0 = RT(yverts[k])
        if opts.norm
            xv[k] = xv0 * RT(beamr)
            yv[k] = yv0 * RT(beamr)
        else
            xv[k] = xv0
            yv[k] = yv0
        end
    end

    ny, nx = size(wf.field)
    size(image) == (ny, nx) || throw(ArgumentError("output size must match wavefront"))
    dx = RT(wf.sampling_m)
    cx = nx ÷ 2
    cy = ny ÷ 2
    xv_dev = backend_adapt(image, xv)
    yv_dev = backend_adapt(image, yv)

    return ka_irregular_polygon_mask!(
        image,
        xv_dev,
        yv_dev,
        cx,
        cy,
        dx;
        dark=opts.dark,
        nsub=antialias_subsampling(),
    )
end

function _prop_irregular_polygon(
    image::AbstractMatrix,
    wf::WaveFront,
    xverts::AbstractVector,
    yverts::AbstractVector,
    opts::IrregularPolygonOptions,
)
    return _prop_irregular_polygon(geometry_exec_style(typeof(image), size(image, 1), size(image, 2)), image, wf, xverts, yverts, opts)
end

function _prop_irregular_polygon(
    wf::WaveFront,
    xverts::AbstractVector,
    yverts::AbstractVector,
    opts::IrregularPolygonOptions,
)
    RT = real(eltype(wf.field))
    image = similar(wf.field, RT, size(wf.field)...)
    return _prop_irregular_polygon(image, wf, xverts, yverts, opts)
end

"""
    prop_irregular_polygon!(image, wf, xverts, yverts; kwargs...)

Write an antialiased convex or concave polygon mask into `image`.

# Arguments
- `image`: destination mask array
- `wf`: wavefront defining the grid and sampling
- `xverts`, `yverts`: polygon vertices

# Keywords
- `DARK`: return the complementary dark mask
- `NORM`: interpret vertices in normalized beam-radius units
"""
function prop_irregular_polygon!(
    image::AbstractMatrix,
    wf::WaveFront,
    xverts::AbstractVector,
    yverts::AbstractVector;
    kwargs...,
)
    return _prop_irregular_polygon(image, wf, xverts, yverts, IrregularPolygonOptions(kwargs))
end

"""Return antialiased filled convex/concave polygon mask."""
function prop_irregular_polygon(
    wf::WaveFront,
    xverts::AbstractVector,
    yverts::AbstractVector;
    kwargs...,
)
    opts = IrregularPolygonOptions(kwargs)
    return _prop_irregular_polygon(wf, xverts, yverts, opts)
end
