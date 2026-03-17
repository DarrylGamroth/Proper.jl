struct MagnifyOptions
    quick::Bool
    conserve::Bool
    amp_conserve::Bool
end

@inline function MagnifyOptions(kwargs::Base.Iterators.Pairs)
    return MagnifyOptions(
        kw_lookup_bool(kwargs, :QUICK, false),
        kw_lookup_bool(kwargs, :CONSERVE, false),
        kw_lookup_bool(kwargs, :AMP_CONSERVE, false),
    )
end

@inline function _apply_magnify_conserve!(out::AbstractMatrix, image_in::AbstractMatrix, mag::Real, opts::MagnifyOptions)
    if opts.conserve
        if eltype(image_in) <: Complex
            out ./= mag
        else
            out ./= mag^2
        end
    elseif opts.amp_conserve
        out ./= mag
    end

    return out
end

@inline function _prop_magnify!(
    out::AbstractMatrix,
    sty::InterpStyle,
    fillsty::AxisFillExecStyle,
    image_in::AbstractMatrix,
    mag::Real,
    opts::MagnifyOptions,
    ws::InterpWorkspace,
    sws::SamplingWorkspace,
)
    ny_out, nx_out = size(out)
    if opts.quick
        T = eltype(ws.xcoords)
        ny, nx = size(image_in)
        cx_in = T(nx ÷ 2)
        cy_in = T(ny ÷ 2)
        cx_out = T(nx_out ÷ 2)
        cy_out = T(ny_out ÷ 2)

        xcoords, ycoords = ensure_interp_axes!(ws, nx_out, ny_out)
        invmag = inv(T(mag))
        fill_affine_axis!(fillsty, xcoords, cx_out, invmag, cx_in)
        fill_affine_axis!(fillsty, ycoords, cy_out, invmag, cy_in)
        prop_cubic_conv_grid!(out, sty, image_in, xcoords, ycoords)
    else
        T = typeof(float(real(zero(eltype(image_in)))))
        _prop_szoom!(out, image_in, T(mag), sws)
    end

    return _apply_magnify_conserve!(out, image_in, mag, opts)
end

@inline function prop_magnify!(
    out::AbstractMatrix,
    image_in::AbstractMatrix,
    mag0::Real,
    opts::MagnifyOptions,
    ctx::RunContext,
)
    return _prop_magnify!(
        out,
        interp_style(ctx),
        axis_fill_exec_style(ctx.backend),
        image_in,
        float(mag0),
        opts,
        interp_workspace(ctx),
        sampling_workspace(ctx),
    )
end

@inline function prop_magnify!(
    out::AbstractMatrix,
    image_in::AbstractMatrix,
    mag0::Real,
    opts::MagnifyOptions,
)
    Tin = typeof(real(zero(eltype(image_in))))
    T = float(promote_type(typeof(mag0), Tin))
    sws = SamplingWorkspace(typeof(out), T)
    return _prop_magnify!(
        out,
        interp_style(typeof(image_in)),
        axis_fill_exec_style(backend_style(typeof(out))),
        image_in,
        float(mag0),
        opts,
        InterpWorkspace(typeof(out), T),
        sws,
    )
end

"""
    prop_magnify!(out, image_in, mag; kwargs...)
    prop_magnify!(out, image_in, mag, ctx; kwargs...)

Resample `image_in` into `out` using either damped-sinc or cubic
interpolation.

# Arguments
- `out`: destination array
- `image_in`: input array to be magnified
- `mag`: magnification factor; for example, `0.5` shrinks the image by a
  factor of two

# Keywords
- `QUICK`: use cubic interpolation instead of the slower damped-sinc path
- `CONSERVE`: conserve intensity; complex arrays are treated as electric
  fields, real arrays as intensity
- `AMP_CONSERVE`: treat a real-valued image as amplitude rather than intensity

# Notes
- The default path uses `prop_szoom` and supports non-square arrays.
- `QUICK=true` uses cubic interpolation and is usually faster but less exact.
"""
function prop_magnify!(
    out::AbstractMatrix,
    image_in::AbstractMatrix,
    mag0::Real;
    kwargs...,
)
    return prop_magnify!(out, image_in, mag0, MagnifyOptions(kwargs))
end

function prop_magnify!(
    out::AbstractMatrix,
    image_in::AbstractMatrix,
    mag0::Real,
    ctx::RunContext;
    kwargs...,
)
    return prop_magnify!(out, image_in, mag0, MagnifyOptions(kwargs), ctx)
end

"""
    prop_magnify(image_in, mag, size_out=0; kwargs...)
    prop_magnify(image_in, mag, size_out, ctx; kwargs...)

Return a magnified or demagnified copy of an image.

# Arguments
- `image_in`: input array
- `mag`: magnification factor

# Keywords
- `size_out`: output dimension when a square output is desired; if zero, the
  output dimensions default to `fix.(size(image_in) .* mag)` for positive
  magnifications
- `QUICK`, `CONSERVE`, `AMP_CONSERVE`: same meanings as in `prop_magnify!`

# Returns
- A new magnified or demagnified array.

# Notes
- The default PROPER behavior is the damped-sinc resampler.
- When `size_out == 0`, Julia computes output height and width independently.
"""
function prop_magnify(
    image_in::AbstractMatrix,
    mag0::Real,
    size_out0::Integer=0;
    kwargs...,
)
    mag = float(mag0)
    opts = MagnifyOptions(kwargs)
    ny, nx = size(image_in)
    if size_out0 > 0
        out_n = Int(size_out0)
        out_n > 0 || throw(ArgumentError("size_out must be positive"))
        out = similar(image_in, out_n, out_n)
        return prop_magnify!(out, image_in, mag, opts)
    end

    if opts.quick
        noy = floor(Int, ny * mag)
        nox = floor(Int, nx * mag)
        noy > 0 || throw(ArgumentError("output height must be positive"))
        nox > 0 || throw(ArgumentError("output width must be positive"))
        out = similar(image_in, noy, nox)
        return prop_magnify!(out, image_in, mag, opts)
    end

    noy = floor(Int, ny * mag)
    nox = floor(Int, nx * mag)
    noy > 0 || throw(ArgumentError("output height must be positive"))
    nox > 0 || throw(ArgumentError("output width must be positive"))
    out = similar(image_in, noy, nox)
    return prop_magnify!(out, image_in, mag, opts)
end

function prop_magnify(
    image_in::AbstractMatrix,
    mag0::Real,
    size_out0::Integer,
    ctx::RunContext;
    kwargs...,
)
    mag = float(mag0)
    opts = MagnifyOptions(kwargs)
    ny, nx = size(image_in)
    if size_out0 > 0
        out_n = Int(size_out0)
        out_n > 0 || throw(ArgumentError("size_out must be positive"))
        out = similar(image_in, out_n, out_n)
        return prop_magnify!(out, image_in, mag, opts, ctx)
    end

    if opts.quick
        noy = floor(Int, ny * mag)
        nox = floor(Int, nx * mag)
        noy > 0 || throw(ArgumentError("output height must be positive"))
        nox > 0 || throw(ArgumentError("output width must be positive"))
        out = similar(image_in, noy, nox)
        return prop_magnify!(out, image_in, mag, opts, ctx)
    end

    noy = floor(Int, ny * mag)
    nox = floor(Int, nx * mag)
    noy > 0 || throw(ArgumentError("output height must be positive"))
    nox > 0 || throw(ArgumentError("output width must be positive"))
    out = similar(image_in, noy, nox)
    return prop_magnify!(out, image_in, mag, opts, ctx)
end
