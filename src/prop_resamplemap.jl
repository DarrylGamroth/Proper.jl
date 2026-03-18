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

Resample `dmap` onto the current wavefront grid and write the result to `out`.

# Arguments
- `out`: destination map with the same size as `wf.field`
- `wf`: wavefront whose grid size and sampling define the output grid
- `dmap`: aberration or amplitude map to be resampled
- `pixscale`: spacing of `dmap` in meters per pixel
- `xc`, `yc`: pixel coordinates of the input map center; `(0, 0)` is the
  center of the first pixel in the upstream convention

# Keywords
- `xshift`, `yshift`: physical shift of the map in meters before resampling

# Notes
- The output dimensions always match `size(wf.field)`.
- The interpolation path follows the accepted executable-baseline cubic
  resampling contract used by parity tests.
- This is the performance-facing entry point for repeated GPU use when paired
  with an explicit `ctx`. The no-`ctx` overload is a convenience wrapper and
  may allocate fresh interpolation state.
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

# Arguments
- `wf`: wavefront defining the target grid
- `dmap`: input map
- `pixscale`: sampling of `dmap` in meters per pixel
- `xc`, `yc`: input-map center in pixel coordinates

# Keywords
- `xshift`, `yshift`: physical map shift in meters

# Returns
- A new array with the same size as `wf.field`.

# Notes
- This is the allocating convenience wrapper around `prop_resamplemap!`.
- For repeated CPU/GPU work, prefer `prop_resamplemap!(out, ..., ctx)` so the
  interpolation workspace is reused.
- The behavior matches the executable PROPER baseline used by the parity
  harnesses.
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
