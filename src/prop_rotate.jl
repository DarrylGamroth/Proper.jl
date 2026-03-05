@enum RotateMethod ROTATE_LINEAR ROTATE_CUBIC

struct RotateOptions{T<:AbstractFloat}
    method::RotateMethod
    cx::T
    cy::T
    sx::T
    sy::T
end

@inline function _rotate_method(kwargs::Base.Iterators.Pairs)::RotateMethod
    m = kw_lookup_string(kwargs, :METH, nothing)
    if m !== nothing
        return lowercase(m) == "linear" ? ROTATE_LINEAR : ROTATE_CUBIC
    end
    if kw_lookup_present(kwargs, :CUBIC)
        return ROTATE_CUBIC
    end
    return ROTATE_CUBIC
end

@inline function RotateOptions(old_image::AbstractMatrix, kwargs::Base.Iterators.Pairs)
    ny, nx = size(old_image)
    method = _rotate_method(kwargs)
    T = float(promote_type(real(eltype(old_image)), Float64))
    cxv = T(kw_lookup_float(kwargs, :XC, nx ÷ 2 + 1))
    cyv = T(kw_lookup_float(kwargs, :YC, ny ÷ 2 + 1))
    sxv = T(kw_lookup_float(kwargs, :XSHIFT, 0.0))
    syv = T(kw_lookup_float(kwargs, :YSHIFT, 0.0))
    return RotateOptions{T}(method, cxv, cyv, sxv, syv)
end

@inline function _prop_rotate_linear!(
    out::AbstractMatrix,
    old_image::AbstractMatrix,
    c::Real,
    s::Real,
    opts::RotateOptions,
)
    ny, nx = size(old_image)
    @inbounds for j in 1:nx
        x = j - opts.cx - opts.sx
        for i in 1:ny
            y = i - opts.cy - opts.sy
            xr = c * x - s * y + opts.cx
            yr = s * x + c * y + opts.cy
            out[i, j] = bilinear_sample(old_image, yr, xr)
        end
    end
    return out
end

@inline function _prop_rotate_cubic!(
    sty::InterpStyle,
    out::AbstractMatrix,
    old_image::AbstractMatrix,
    c::Real,
    s::Real,
    opts::RotateOptions,
)
    ny, nx = size(old_image)
    @inbounds for j in 1:nx
        x = j - opts.cx - opts.sx
        for i in 1:ny
            y = i - opts.cy - opts.sy
            xr = c * x - s * y + opts.cx
            yr = s * x + c * y + opts.cy
            out[i, j] = prop_cubic_conv(sty, old_image, yr, xr)
        end
    end
    return out
end

@inline function _prop_rotate!(
    sty::InterpStyle,
    out::AbstractMatrix,
    old_image::AbstractMatrix,
    theta::Real,
    opts::RotateOptions,
)
    size(out) == size(old_image) || throw(ArgumentError("output size must match input image"))
    if !(old_image isa StridedMatrix && out isa StridedMatrix)
        host_out = Matrix{eltype(out)}(undef, size(out)...)
        _prop_rotate!(sty, host_out, Matrix(old_image), theta, opts)
        copyto!(out, host_out)
        return out
    end

    ang = deg2rad(-float(theta))
    c = cos(ang)
    s = sin(ang)
    return opts.method === ROTATE_LINEAR ?
        _prop_rotate_linear!(out, old_image, c, s, opts) :
        _prop_rotate_cubic!(sty, out, old_image, c, s, opts)
end

@inline function prop_rotate!(out::AbstractMatrix, old_image::AbstractMatrix, theta::Real, opts::RotateOptions, ctx::RunContext)
    return _prop_rotate!(interp_style(ctx), out, old_image, theta, opts)
end

@inline function prop_rotate!(out::AbstractMatrix, old_image::AbstractMatrix, theta::Real, opts::RotateOptions)
    return _prop_rotate!(interp_style(typeof(old_image)), out, old_image, theta, opts)
end

@inline function prop_rotate!(out::AbstractMatrix, old_image::AbstractMatrix, theta::Real; kwargs...)
    return prop_rotate!(out, old_image, theta, RotateOptions(old_image, kwargs))
end

@inline function prop_rotate!(out::AbstractMatrix, old_image::AbstractMatrix, theta::Real, ctx::RunContext; kwargs...)
    return prop_rotate!(out, old_image, theta, RotateOptions(old_image, kwargs), ctx)
end

@inline function prop_rotate(old_image::AbstractMatrix, theta::Real, opts::RotateOptions, ctx::RunContext)
    out = similar(old_image)
    return prop_rotate!(out, old_image, theta, opts, ctx)
end

@inline function prop_rotate(old_image::AbstractMatrix, theta::Real, opts::RotateOptions)
    out = similar(old_image)
    return prop_rotate!(out, old_image, theta, opts)
end

"""Rotate image by theta degrees around center (cubic by default)."""
function prop_rotate(
    old_image::AbstractMatrix,
    theta::Real;
    kwargs...,
)
    opts = RotateOptions(old_image, kwargs)
    return prop_rotate(old_image, theta, opts)
end

function prop_rotate(old_image::AbstractMatrix, theta::Real, ctx::RunContext; kwargs...)
    opts = RotateOptions(old_image, kwargs)
    return prop_rotate(old_image, theta, opts, ctx)
end
