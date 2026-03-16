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
        prop_szoom!(out, image_in, mag)
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
    return _prop_magnify!(
        out,
        interp_style(typeof(image_in)),
        axis_fill_exec_style(backend_style(typeof(out))),
        image_in,
        float(mag0),
        opts,
        InterpWorkspace(typeof(out), T),
    )
end

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

"""Magnify image using upstream PROPER behavior (`QUICK` cubic, default damped-sinc)."""
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
