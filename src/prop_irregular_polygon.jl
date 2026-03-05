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

"""Return antialiased filled convex/concave polygon mask."""
function prop_irregular_polygon(wf::WaveFront, xverts::AbstractVector, yverts::AbstractVector; kwargs...)
    length(xverts) == length(yverts) || throw(ArgumentError("vertex arrays must have same length"))

    dark = switch_set(:DARK; kwargs...)
    norm = switch_set(:NORM; kwargs...)
    beamr = prop_get_beamradius(wf)

    xv = norm ? float.(xverts) .* beamr : float.(xverts)
    yv = norm ? float.(yverts) .* beamr : float.(yverts)

    ny, nx = size(wf.field)
    dx = wf.sampling_m
    x = coordinate_axis(nx, dx)
    y = coordinate_axis(ny, dx)

    subs = 3
    offs = ((collect(1:subs) .- (subs + 1) / 2) ./ subs) .* dx

    image = zeros(Float64, ny, nx)
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

    return dark ? (1 .- image) : image
end
