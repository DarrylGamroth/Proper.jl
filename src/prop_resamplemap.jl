struct ResampleMapOptions{T<:AbstractFloat}
    pixscale::T
    xc::T
    yc::T
    xshift::T
    yshift::T
end

@inline function ResampleMapOptions(
    wf::WaveFront,
    pixscale::Real,
    xc::Real,
    yc::Real,
    xshift::Real=0.0,
    yshift::Real=0.0,
)
    T = float(promote_type(typeof(pixscale), typeof(xc), typeof(yc), typeof(xshift), typeof(yshift), typeof(wf.sampling_m)))
    return ResampleMapOptions{T}(T(pixscale), T(xc), T(yc), T(xshift), T(yshift))
end

"""Resample map to wavefront grid with cubic interpolation and physical shifts."""
function prop_resamplemap!(
    out::AbstractMatrix,
    wf::WaveFront,
    dmap::AbstractMatrix,
    opts::ResampleMapOptions,
)
    ny, nx = size(wf.field)
    size(out) == (ny, nx) || throw(ArgumentError("output size must match wavefront"))

    T = typeof(opts.pixscale)
    scale = T(wf.sampling_m) / opts.pixscale
    xoff = opts.xc - opts.xshift / opts.pixscale
    yoff = opts.yc - opts.yshift / opts.pixscale

    xcoords = Vector{T}(undef, nx)
    ycoords = Vector{T}(undef, ny)

    @inbounds for j in 1:nx
        xcoords[j] = (T(j - 1 - (nx ÷ 2)) * scale) + xoff
    end
    @inbounds for i in 1:ny
        ycoords[i] = (T(i - 1 - (ny ÷ 2)) * scale) + yoff
    end

    sampled = prop_cubic_conv(dmap, xcoords, ycoords; grid=true)
    copyto!(out, sampled)
    return out
end

@inline function prop_resamplemap!(
    out::AbstractMatrix,
    wf::WaveFront,
    dmap::AbstractMatrix,
    pixscale::Real,
    xc::Real,
    yc::Real,
    xshift::Real=0.0,
    yshift::Real=0.0,
)
    return prop_resamplemap!(out, wf, dmap, ResampleMapOptions(wf, pixscale, xc, yc, xshift, yshift))
end

function prop_resamplemap(wf::WaveFront, dmap::AbstractMatrix, pixscale::Real, xc::Real, yc::Real, xshift::Real=0.0, yshift::Real=0.0)
    ny, nx = size(wf.field)
    Tout = float(promote_type(real(eltype(dmap)), typeof(wf.sampling_m), typeof(pixscale), typeof(xc), typeof(yc), typeof(xshift), typeof(yshift)))
    out = similar(dmap, Tout, ny, nx)
    return prop_resamplemap!(out, wf, dmap, pixscale, xc, yc, xshift, yshift)
end
