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

@inline function _prop_magnify(sty::InterpStyle, image_in::AbstractMatrix, mag::Real, out_n::Int, opts::MagnifyOptions)
    out = if opts.quick
        Tin = typeof(real(zero(eltype(image_in))))
        T = float(promote_type(typeof(mag), Tin))
        ny, nx = size(image_in)
        cx_in = T(nx ÷ 2)
        cy_in = T(ny ÷ 2)
        cx_out = T(out_n ÷ 2)
        cy_out = T(out_n ÷ 2)

        xcoords = Vector{T}(undef, out_n)
        ycoords = Vector{T}(undef, out_n)
        invmag = inv(T(mag))

        @inbounds for j in 1:out_n
            xcoords[j] = (T(j - 1) - cx_out) * invmag + cx_in
        end
        @inbounds for i in 1:out_n
            ycoords[i] = (T(i - 1) - cy_out) * invmag + cy_in
        end
        prop_cubic_conv(sty, image_in, xcoords, ycoords; grid=true)
    else
        prop_szoom(image_in, mag, out_n)
    end

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

@inline function _prop_magnify(image_in::AbstractMatrix, mag::Real, out_n::Int, opts::MagnifyOptions)
    return _prop_magnify(interp_style(typeof(image_in)), image_in, mag, out_n, opts)
end

@inline function _prop_magnify(image_in::AbstractMatrix, mag::Real, out_n::Int, opts::MagnifyOptions, ctx::RunContext)
    return _prop_magnify(interp_style(ctx), image_in, mag, out_n, opts)
end

"""Magnify image using upstream PROPER behavior (`QUICK` cubic, default damped-sinc)."""
function prop_magnify(
    image_in::AbstractMatrix,
    mag0::Real,
    size_out0::Integer=0;
    kwargs...,
)
    mag = float(mag0)
    ny, _ = size(image_in)
    out_n = size_out0 > 0 ? Int(size_out0) : round(Int, ny * mag)
    out_n > 0 || throw(ArgumentError("size_out must be positive"))
    opts = MagnifyOptions(kwargs)
    return _prop_magnify(image_in, mag, out_n, opts)
end

function prop_magnify(
    image_in::AbstractMatrix,
    mag0::Real,
    size_out0::Integer,
    ctx::RunContext;
    kwargs...,
)
    mag = float(mag0)
    ny, _ = size(image_in)
    out_n = size_out0 > 0 ? Int(size_out0) : round(Int, ny * mag)
    out_n > 0 || throw(ArgumentError("size_out must be positive"))
    opts = MagnifyOptions(kwargs)
    return _prop_magnify(image_in, mag, out_n, opts, ctx)
end
