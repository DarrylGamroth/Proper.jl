@inline function _point_in_poly(x::Real, y::Real, xv::AbstractVector, yv::AbstractVector)
    inside = false
    n = length(xv)
    j = n
    @inbounds for i in 1:n
        yi = yv[i]
        yj = yv[j]
        if ((yi > y) != (yj > y))
            xcross = (xv[j] - xv[i]) * (y - yi) / (yj - yi + eps()) + xv[i]
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
    image::AbstractMatrix,
    wf::WaveFront,
    xverts::AbstractVector,
    yverts::AbstractVector,
    opts::IrregularPolygonOptions,
)
    length(xverts) == length(yverts) || throw(ArgumentError("vertex arrays must have same length"))

    beamr = prop_get_beamradius(wf)
    xv = opts.norm ? float.(xverts) .* beamr : float.(xverts)
    yv = opts.norm ? float.(yverts) .* beamr : float.(yverts)

    ny, nx = size(wf.field)
    size(image) == (ny, nx) || throw(ArgumentError("output size must match wavefront"))
    dx = wf.sampling_m
    x = coordinate_axis(nx, dx)
    y = coordinate_axis(ny, dx)

    subs = antialias_subsampling()
    offs = ((collect(1:subs) .- (subs + 1) / 2) ./ subs) .* dx

    fill!(image, zero(eltype(image)))
    @inbounds for j in 1:nx
        x0 = x[j]
        for i in 1:ny
            y0 = y[i]
            cnt = 0
            for oy in offs, ox in offs
                cnt += _point_in_poly(x0 + ox, y0 + oy, xv, yv) ? 1 : 0
            end
            image[i, j] = cnt / (subs * subs)
        end
    end

    if opts.dark
        image .= 1 .- image
    end
    return image
end

function _prop_irregular_polygon(
    wf::WaveFront,
    xverts::AbstractVector,
    yverts::AbstractVector,
    opts::IrregularPolygonOptions,
)
    RT = real(eltype(wf.field))
    image = zeros(RT, size(wf.field)...)
    return _prop_irregular_polygon(image, wf, xverts, yverts, opts)
end

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
