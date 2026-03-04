"""Generate binary rectangular mask centered at (xc, yc) in meters."""
function prop_rectangle(wf::WaveFront, xsize::Real, ysize::Real, xc::Real=0.0, yc::Real=0.0; kwargs...)
    ny, nx = size(wf.field)
    dx = wf.sampling_m
    x = coordinate_axis(nx, dx) .- float(xc)
    y = coordinate_axis(ny, dx) .- float(yc)
    hx = float(xsize) / 2
    hy = float(ysize) / 2
    mask = similar(real.(wf.field), ny, nx)
    @inbounds for j in 1:nx
        xin = abs(x[j]) <= hx
        for i in 1:ny
            mask[i, j] = (xin && abs(y[i]) <= hy) ? one(eltype(mask)) : zero(eltype(mask))
        end
    end
    return mask
end
