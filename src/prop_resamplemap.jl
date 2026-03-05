"""Resample map to wavefront grid with cubic interpolation and physical shifts."""
function prop_resamplemap(wf::WaveFront, dmap::AbstractMatrix, pixscale::Real, xc::Real, yc::Real, xshift::Real=0.0, yshift::Real=0.0)
    ny, nx = size(wf.field)
    T = float(promote_type(typeof(pixscale), typeof(xc), typeof(yc), typeof(xshift), typeof(yshift), typeof(wf.sampling_m)))

    scale = T(wf.sampling_m) / T(pixscale)
    xoff = T(xc) - T(xshift) / T(pixscale)
    yoff = T(yc) - T(yshift) / T(pixscale)

    xcoords = Vector{T}(undef, nx)
    ycoords = Vector{T}(undef, ny)

    @inbounds for j in 1:nx
        xcoords[j] = (T(j - 1 - (nx ÷ 2)) * scale) + xoff
    end
    @inbounds for i in 1:ny
        ycoords[i] = (T(i - 1 - (ny ÷ 2)) * scale) + yoff
    end

    return prop_cubic_conv(dmap, xcoords, ycoords; grid=true)
end
