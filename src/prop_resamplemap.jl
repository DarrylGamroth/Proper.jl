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
function _prop_resamplemap!(
    sty::InterpStyle,
    out::AbstractMatrix,
    wf::WaveFront,
    dmap::AbstractMatrix,
    opts::ResampleMapOptions,
    ws::InterpWorkspace,
)
    ny, nx = size(wf.field)
    size(out) == (ny, nx) || throw(ArgumentError("output size must match wavefront"))

    xcoords, ycoords = ensure_interp_axes!(ws, nx, ny)
    T = eltype(xcoords)
    scale = T(wf.sampling_m) / T(opts.pixscale)
    xoff = T(opts.xc) - T(opts.xshift) / T(opts.pixscale)
    yoff = T(opts.yc) - T(opts.yshift) / T(opts.pixscale)

    @inbounds for j in 1:nx
        xcoords[j] = (T(j - 1 - (nx ÷ 2)) * scale) + xoff
    end
    @inbounds for i in 1:ny
        ycoords[i] = (T(i - 1 - (ny ÷ 2)) * scale) + yoff
    end

    return prop_cubic_conv_grid!(out, sty, dmap, xcoords, ycoords)
end

function prop_resamplemap!(
    out::AbstractMatrix,
    wf::WaveFront,
    dmap::AbstractMatrix,
    opts::ResampleMapOptions,
    ctx::RunContext,
)
    return _prop_resamplemap!(interp_style(ctx), out, wf, dmap, opts, interp_workspace(ctx))
end

function prop_resamplemap!(
    out::AbstractMatrix,
    wf::WaveFront,
    dmap::AbstractMatrix,
    opts::ResampleMapOptions,
)
    return _prop_resamplemap!(
        interp_style(typeof(out)),
        out,
        wf,
        dmap,
        opts,
        InterpWorkspace(float(typeof(opts.pixscale))),
    )
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
    return prop_resamplemap!(out, wf, dmap, ResampleMapOptions(wf, pixscale, xc, yc, xshift, yshift), RunContext(typeof(out)))
end

@inline function prop_resamplemap!(
    out::AbstractMatrix,
    wf::WaveFront,
    dmap::AbstractMatrix,
    pixscale::Real,
    xc::Real,
    yc::Real,
    ctx::RunContext,
    xshift::Real=0.0,
    yshift::Real=0.0,
)
    return prop_resamplemap!(out, wf, dmap, ResampleMapOptions(wf, pixscale, xc, yc, xshift, yshift), ctx)
end

function prop_resamplemap(wf::WaveFront, dmap::AbstractMatrix, pixscale::Real, xc::Real, yc::Real, xshift::Real=0.0, yshift::Real=0.0)
    ny, nx = size(wf.field)
    Tout = float(promote_type(real(eltype(dmap)), typeof(wf.sampling_m), typeof(pixscale), typeof(xc), typeof(yc), typeof(xshift), typeof(yshift)))
    out = similar(dmap, Tout, ny, nx)
    return prop_resamplemap!(out, wf, dmap, pixscale, xc, yc, RunContext(typeof(out)), xshift, yshift)
end

function prop_resamplemap(wf::WaveFront, dmap::AbstractMatrix, pixscale::Real, xc::Real, yc::Real, ctx::RunContext, xshift::Real=0.0, yshift::Real=0.0)
    ny, nx = size(wf.field)
    Tout = float(promote_type(real(eltype(dmap)), typeof(wf.sampling_m), typeof(pixscale), typeof(xc), typeof(yc), typeof(xshift), typeof(yshift)))
    out = similar(dmap, Tout, ny, nx)
    return prop_resamplemap!(out, wf, dmap, pixscale, xc, yc, ctx, xshift, yshift)
end
