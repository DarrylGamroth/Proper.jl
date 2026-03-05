struct PolygonOptions{T<:AbstractFloat}
    norm::Bool
    dark::Bool
    rotation::T
end

@inline function PolygonOptions(kwargs::Base.Iterators.Pairs)
    return PolygonOptions{Float64}(
        kw_lookup_bool(kwargs, :NORM, false),
        kw_lookup_bool(kwargs, :DARK, false),
        kw_lookup_float(kwargs, :ROTATION, 0.0),
    )
end

function _prop_polygon(wf::WaveFront, nsides::Integer, radius::Real, xc::Real, yc::Real, opts::PolygonOptions)
    nsides >= 3 || throw(ArgumentError("nsides must be >= 3"))

    beamr = prop_get_beamradius(wf)
    rad = opts.norm ? float(radius) * beamr : float(radius)
    x0 = opts.norm ? float(xc) * beamr : float(xc)
    y0 = opts.norm ? float(yc) * beamr : float(yc)
    θ0 = deg2rad(opts.rotation)

    t = -((collect(0:(nsides - 1)) ./ nsides) .* 2pi) .+ θ0
    xv = x0 .+ rad .* cos.(t)
    yv = y0 .+ rad .* sin.(t)

    return _prop_irregular_polygon(wf, xv, yv, IrregularPolygonOptions(opts.dark, false))
end

function _prop_polygon!(
    image::AbstractMatrix,
    wf::WaveFront,
    nsides::Integer,
    radius::Real,
    xc::Real,
    yc::Real,
    opts::PolygonOptions,
)
    nsides >= 3 || throw(ArgumentError("nsides must be >= 3"))

    beamr = prop_get_beamradius(wf)
    rad = opts.norm ? float(radius) * beamr : float(radius)
    x0 = opts.norm ? float(xc) * beamr : float(xc)
    y0 = opts.norm ? float(yc) * beamr : float(yc)
    θ0 = deg2rad(opts.rotation)

    t = -((collect(0:(nsides - 1)) ./ nsides) .* 2pi) .+ θ0
    xv = x0 .+ rad .* cos.(t)
    yv = y0 .+ rad .* sin.(t)

    return _prop_irregular_polygon(image, wf, xv, yv, IrregularPolygonOptions(opts.dark, false))
end

function prop_polygon!(
    image::AbstractMatrix,
    wf::WaveFront,
    nsides::Integer,
    radius::Real,
    xc::Real=0.0,
    yc::Real=0.0;
    kwargs...,
)
    opts = PolygonOptions(kwargs)
    return _prop_polygon!(image, wf, nsides, radius, xc, yc, opts)
end

"""Return antialiased regular polygon mask."""
function prop_polygon(
    wf::WaveFront,
    nsides::Integer,
    radius::Real,
    xc::Real=0.0,
    yc::Real=0.0;
    kwargs...,
)
    opts = PolygonOptions(kwargs)
    return _prop_polygon(wf, nsides, radius, xc, yc, opts)
end
