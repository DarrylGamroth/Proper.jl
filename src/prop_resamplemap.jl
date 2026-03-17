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

"""Internal implementation for map resampling onto the current wavefront grid."""
function _prop_resamplemap!(
    sty::InterpStyle,
    fillsty::AxisFillExecStyle,
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

    fill_affine_axis!(fillsty, xcoords, T(nx ÷ 2), scale, xoff)
    fill_affine_axis!(fillsty, ycoords, T(ny ÷ 2), scale, yoff)

    return prop_cubic_conv_grid!(out, sty, dmap, xcoords, ycoords)
end

"""
    prop_resamplemap!(out, wf, dmap, opts, ctx)
    prop_resamplemap!(out, wf, dmap, pixscale, xc, yc, xshift=0, yshift=0)
    prop_resamplemap!(out, wf, dmap, pixscale, xc, yc, ctx, xshift=0, yshift=0)

Resample an input map onto the current wavefront grid using cubic
interpolation.

The output size and sampling match `wf.field`. `pixscale` is the input map
sampling in meters per pixel. `xc` and `yc` specify the map center in pixel
coordinates. `xshift` and `yshift` are physical shifts in meters applied before
resampling.
"""
function prop_resamplemap!(
    out::AbstractMatrix,
    wf::WaveFront,
    dmap::AbstractMatrix,
    opts::ResampleMapOptions,
    ctx::RunContext,
)
    return _prop_resamplemap!(
        interp_style(ctx),
        axis_fill_exec_style(ctx.backend),
        out,
        wf,
        dmap,
        opts,
        interp_workspace(ctx),
    )
end

function prop_resamplemap!(
    out::AbstractMatrix,
    wf::WaveFront,
    dmap::AbstractMatrix,
    opts::ResampleMapOptions,
)
    return _prop_resamplemap!(
        interp_style(typeof(out)),
        axis_fill_exec_style(backend_style(typeof(out))),
        out,
        wf,
        dmap,
        opts,
        InterpWorkspace(typeof(out), float(typeof(opts.pixscale))),
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

"""
    prop_resamplemap(wf, dmap, pixscale, xc, yc, xshift=0, yshift=0)
    prop_resamplemap(wf, dmap, pixscale, xc, yc, ctx, xshift=0, yshift=0)

Return a new map resampled to the current wavefront grid.

This matches the executable PROPER baseline used by the parity harnesses. The
result has the same dimensions as `wf.field` and the wavefront's current
sampling in meters per pixel.
"""
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
