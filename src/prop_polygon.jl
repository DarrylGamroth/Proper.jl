"""Return antialiased regular polygon mask."""
function prop_polygon(wf::WaveFront, nsides::Integer, radius::Real, xc::Real=0.0, yc::Real=0.0; kwargs...)
    nsides >= 3 || throw(ArgumentError("nsides must be >= 3"))

    norm = switch_set(:NORM; kwargs...)
    beamr = prop_get_beamradius(wf)
    rad = norm ? float(radius) * beamr : float(radius)
    x0 = norm ? float(xc) * beamr : float(xc)
    y0 = norm ? float(yc) * beamr : float(yc)

    rotation = haskey(kwargs, :ROTATION) ? float(kwargs[:ROTATION]) : haskey(kwargs, :rotation) ? float(kwargs[:rotation]) : 0.0
    θ0 = deg2rad(rotation)

    t = -((collect(0:(nsides - 1)) ./ nsides) .* 2pi) .+ θ0
    xv = x0 .+ rad .* cos.(t)
    yv = y0 .+ rad .* sin.(t)

    return prop_irregular_polygon(wf, xv, yv; kwargs...)
end
