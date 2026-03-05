"""Generate binary rectangular mask centered at (xc, yc) in meters."""
function prop_rectangle(wf::WaveFront, xsize::Real, ysize::Real, xc::Real=0.0, yc::Real=0.0; kwargs...)
    norm = switch_set(:NORM; kwargs...)
    beamr = prop_get_beamradius(wf)
    scale = norm ? beamr : 1.0

    ny, nx = size(wf.field)
    dx = wf.sampling_m
    x = coordinate_axis(nx, dx) .- float(xc) * scale
    y = coordinate_axis(ny, dx) .- float(yc) * scale
    hx = float(xsize) * scale / 2
    hy = float(ysize) * scale / 2
    RT = real(eltype(wf.field))
    mask = similar(wf.field, RT, ny, nx)
    @inbounds for j in 1:nx
        xin = abs(x[j]) <= hx
        for i in 1:ny
            mask[i, j] = (xin && abs(y[i]) <= hy) ? one(eltype(mask)) : zero(eltype(mask))
        end
    end
    return mask
end
