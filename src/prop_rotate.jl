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
            out[i, j] = prop_cubic_conv(old_image, yr, xr)
        end
    end
    return out
end

@inline function prop_rotate(old_image::AbstractMatrix, theta::Real, opts::RotateOptions)
    out = similar(old_image)
    ang = deg2rad(-float(theta))
    c = cos(ang)
    s = sin(ang)
    return opts.method === ROTATE_LINEAR ?
        _prop_rotate_linear!(out, old_image, c, s, opts) :
        _prop_rotate_cubic!(out, old_image, c, s, opts)
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
