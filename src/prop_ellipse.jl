"""Generate binary elliptical mask centered at (xc, yc) in meters."""
function prop_ellipse(wf::WaveFront, xradius::Real, yradius::Real, xc::Real=0.0, yc::Real=0.0; kwargs...)
    ny, nx = size(wf.field)
    dx = wf.sampling_m
    x = coordinate_axis(nx, dx) .- float(xc)
    y = coordinate_axis(ny, dx) .- float(yc)
    xr = float(xradius)
    yr = float(yradius)
    mask = similar(real.(wf.field), ny, nx)
    @inbounds for j in 1:nx
        xx = (x[j] / xr)^2
        for i in 1:ny
            mask[i, j] = (xx + (y[i] / yr)^2) <= 1 ? one(eltype(mask)) : zero(eltype(mask))
        end
    end
    return mask
end
