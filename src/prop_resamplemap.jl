"""Resample map to wavefront grid with map-center and physical shifts."""
function prop_resamplemap(wf::WaveFront, dmap::AbstractMatrix, pixscale::Real, xc::Real, yc::Real, xshift::Real=0.0, yshift::Real=0.0)
    ny, nx = size(wf.field)
    out = similar(dmap, ny, nx)

    map_pix_per_m = inv(float(pixscale))
    dx_out = wf.sampling_m

    cx_out = (nx + 1) / 2
    cy_out = (ny + 1) / 2

    @inbounds for j in 1:nx
        x_m = (j - cx_out) * dx_out + float(xshift)
        x_map = float(xc) + x_m * map_pix_per_m
        for i in 1:ny
            y_m = (i - cy_out) * dx_out + float(yshift)
            y_map = float(yc) + y_m * map_pix_per_m
            out[i, j] = bilinear_sample(dmap, y_map, x_map)
        end
    end

    return out
end
