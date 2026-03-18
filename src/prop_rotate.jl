@enum RotateMethod ROTATE_LINEAR ROTATE_CUBIC

struct RotateOptions{T<:AbstractFloat}
    method::RotateMethod
    cx::T
    cy::T
    sx::T
    sy::T
    missing::T
end

@inline function _rotate_method(kwargs::Base.Iterators.Pairs)::RotateMethod
    m = kw_lookup_string(kwargs, :METH, nothing)
    if m !== nothing
        return lowercase(m) == "linear" ? ROTATE_LINEAR : ROTATE_CUBIC
    end
    if kw_lookup_present(kwargs, :CUBIC)
        return ROTATE_CUBIC
    end
    return ROTATE_LINEAR
end

@inline function RotateOptions(old_image::AbstractMatrix, kwargs::Base.Iterators.Pairs)
    ny, nx = size(old_image)
    method = _rotate_method(kwargs)
    T = float(promote_type(real(eltype(old_image)), Float64))
    cxv = T(kw_lookup_float(kwargs, :XC, nx ÷ 2 + 1))
    cyv = T(kw_lookup_float(kwargs, :YC, ny ÷ 2 + 1))
    sxv = T(kw_lookup_float(kwargs, :XSHIFT, 0.0))
    syv = T(kw_lookup_float(kwargs, :YSHIFT, 0.0))
    missingv = T(kw_lookup_float(kwargs, :MISSING, kw_lookup_float(kwargs, :EXTR, 0.0)))
    return RotateOptions{T}(method, cxv, cyv, sxv, syv, missingv)
end

@inline function RotateOptions(old_image::AbstractMatrix)
    ny, nx = size(old_image)
    T = float(promote_type(real(eltype(old_image)), Float64))
    return RotateOptions{T}(ROTATE_LINEAR, T(nx ÷ 2 + 1), T(ny ÷ 2 + 1), zero(T), zero(T), zero(T))
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
            out[i, j] = bilinear_sample_extrap(old_image, yr, xr, opts.missing)
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
            if 1 <= xr <= nx && 1 <= yr <= ny
                # `libcconv` follows the upstream 0-based cubic-convolution kernel,
                # so rotate coordinates must be shifted out of Julia/MATLAB 1-based space.
                out[i, j] = prop_cubic_conv(sty, old_image, xr - 1, yr - 1)
            else
                out[i, j] = opts.missing
            end
        end
    end
    return out
end

abstract type RotateExecStyle end
struct RotateLoopExecStyle <: RotateExecStyle end
struct RotateKAExecStyle <: RotateExecStyle end
struct RotateHostExecStyle <: RotateExecStyle end
struct RotateUnsupportedExecStyle <: RotateExecStyle end

@inline rotate_ka_support(::Val{ROTATE_LINEAR}, ::InterpStyle) = Val(true)
@inline rotate_ka_support(::Val{ROTATE_CUBIC}, ::CubicInterpStyle) = Val(true)
@inline rotate_ka_support(::Val{ROTATE_CUBIC}, ::InterpStyle) = Val(false)

@inline rotate_exec_style(
    ::StridedLayout,
    ::StridedLayout,
    ::CPUBackend,
    ::CPUBackend,
    ::Val{false},
    ::Val,
) = RotateLoopExecStyle()

@inline rotate_exec_style(
    ::StridedLayout,
    ::StridedLayout,
    ::CPUBackend,
    ::CPUBackend,
    ::Val{true},
    ::Val{false},
) = RotateLoopExecStyle()

@inline rotate_exec_style(
    ::ArrayLayoutStyle,
    ::ArrayLayoutStyle,
    ::B,
    ::B,
    ::Val{true},
    ::Val{true},
) where {B<:BackendStyle} = RotateKAExecStyle()

@inline rotate_exec_style(
    ::ArrayLayoutStyle,
    ::ArrayLayoutStyle,
    ::CPUBackend,
    ::CPUBackend,
    ::Val,
    ::Val,
) = RotateHostExecStyle()

@inline rotate_exec_style(
    ::ArrayLayoutStyle,
    ::ArrayLayoutStyle,
    ::BackendStyle,
    ::BackendStyle,
    ::Val,
    ::Val,
) = RotateUnsupportedExecStyle()

@inline function _prop_rotate_ka!(
    ::Val{ROTATE_LINEAR},
    out::AbstractMatrix,
    old_image::AbstractMatrix,
    c::Real,
    s::Real,
    opts::RotateOptions,
    sty::InterpStyle,
)
    _ = sty
    return ka_rotate_linear!(out, old_image, c, s, opts.cx, opts.cy, opts.sx, opts.sy, opts.missing)
end

@inline function _prop_rotate_ka!(
    ::Val{ROTATE_CUBIC},
    out::AbstractMatrix,
    old_image::AbstractMatrix,
    c::Real,
    s::Real,
    opts::RotateOptions,
    ::CubicInterpStyle,
)
    return ka_rotate_cubic!(out, old_image, c, s, opts.cx, opts.cy, opts.sx, opts.sy, opts.missing)
end

@inline function _prop_rotate_loop!(
    ::Val{ROTATE_LINEAR},
    sty::InterpStyle,
    out::StridedMatrix,
    old_image::StridedMatrix,
    c::Real,
    s::Real,
    opts::RotateOptions,
)
    _ = sty
    return _prop_rotate_linear!(out, old_image, c, s, opts)
end

@inline function _prop_rotate_loop!(
    ::Val{ROTATE_CUBIC},
    sty::InterpStyle,
    out::StridedMatrix,
    old_image::StridedMatrix,
    c::Real,
    s::Real,
    opts::RotateOptions,
)
    return _prop_rotate_cubic!(sty, out, old_image, c, s, opts)
end

@inline function _prop_rotate_exec!(
    ::RotateKAExecStyle,
    sty::InterpStyle,
    out::AbstractMatrix,
    old_image::AbstractMatrix,
    theta::Real,
    opts::RotateOptions,
)
    ang = deg2rad(-float(theta))
    c = cos(ang)
    s = sin(ang)
    return _prop_rotate_ka!(Val(opts.method), out, old_image, c, s, opts, sty)
end

@inline function _prop_rotate_exec!(
    ::RotateLoopExecStyle,
    sty::InterpStyle,
    out::StridedMatrix,
    old_image::StridedMatrix,
    theta::Real,
    opts::RotateOptions,
)
    ang = deg2rad(-float(theta))
    c = cos(ang)
    s = sin(ang)
    return _prop_rotate_loop!(Val(opts.method), sty, out, old_image, c, s, opts)
end

@inline function _prop_rotate_exec!(
    ::RotateHostExecStyle,
    sty::InterpStyle,
    out::AbstractMatrix,
    old_image::AbstractMatrix,
    theta::Real,
    opts::RotateOptions,
)
    host_out = Matrix{eltype(out)}(undef, size(out)...)
    _prop_rotate!(sty, host_out, Matrix(old_image), theta, opts)
    copyto!(out, host_out)
    return out
end

@inline function _prop_rotate_exec!(
    ::RotateUnsupportedExecStyle,
    sty::InterpStyle,
    out::AbstractMatrix,
    old_image::AbstractMatrix,
    theta::Real,
    opts::RotateOptions,
)
    _ = sty
    _ = theta
    _ = opts
    throw(ArgumentError("prop_rotate has no native implementation for $(typeof(old_image)) -> $(typeof(out)); use matching CPU arrays or a backend-native path"))
end

@inline function _prop_rotate!(
    sty::InterpStyle,
    out::AbstractMatrix,
    old_image::AbstractMatrix,
    theta::Real,
    opts::RotateOptions,
)
    size(out) == size(old_image) || throw(ArgumentError("output size must match input image"))
    ny, nx = size(out)
    sty_exec = rotate_exec_style(
        array_layout_style(typeof(out)),
        array_layout_style(typeof(old_image)),
        backend_style(typeof(out)),
        backend_style(typeof(old_image)),
        Val(ka_rotate_enabled(typeof(out), ny, nx)),
        rotate_ka_support(Val(opts.method), sty),
    )
    return _prop_rotate_exec!(sty_exec, sty, out, old_image, theta, opts)
end

"""
    prop_rotate!(out, old_image, theta, opts, ctx)
    prop_rotate!(out, old_image, theta, opts)
    prop_rotate!(out, old_image, theta; kwargs...)
    prop_rotate!(out, old_image, theta, ctx; kwargs...)

Rotate `old_image` by `theta` degrees counter-clockwise and write the result to
`out`.

# Arguments
- `out`: destination array with the same size as `old_image`
- `old_image`: input image
- `theta`: rotation angle in degrees counter-clockwise

# Keywords
- `XC`, `YC`: pixel coordinates of the input-array center; defaults are
  `fix(nx/2)+1` and `fix(ny/2)+1`
- `XSHIFT`, `YSHIFT`: shift of the output image in pixels
- `MISSING` / `EXTR`: value assigned to extrapolated pixels, default `0`
- `METH`: interpolation method; accepted public values are `"linear"` and
  `"cubic"`
- `CUBIC`: compatibility shortcut for `METH="cubic"`

# Notes
- The accepted Julia semantics follow the MATLAB-aligned decision for the
  default path: linear interpolation unless cubic is explicitly requested.
- The cubic path uses the upstream PROPER cubic-convolution kernel rather than
  a general spline/interp2 implementation.
- GPU/backend performance paths require native backend support; unsupported
  backend combinations throw instead of silently materializing on the CPU.
"""
@inline function prop_rotate!(out::AbstractMatrix, old_image::AbstractMatrix, theta::Real, opts::RotateOptions, ctx::RunContext)
    return _prop_rotate!(interp_style(ctx), out, old_image, theta, opts)
end

@inline function prop_rotate!(out::AbstractMatrix, old_image::AbstractMatrix, theta::Real, opts::RotateOptions)
    return _prop_rotate!(interp_style(typeof(old_image)), out, old_image, theta, opts)
end

@inline function prop_rotate!(out::AbstractMatrix, old_image::AbstractMatrix, theta::Real; kwargs...)
    opts = isempty(kwargs) ? RotateOptions(old_image) : RotateOptions(old_image, kwargs)
    return prop_rotate!(out, old_image, theta, opts)
end

@inline function prop_rotate!(out::AbstractMatrix, old_image::AbstractMatrix, theta::Real, ctx::RunContext; kwargs...)
    opts = isempty(kwargs) ? RotateOptions(old_image) : RotateOptions(old_image, kwargs)
    return prop_rotate!(out, old_image, theta, opts, ctx)
end

@inline function prop_rotate(old_image::AbstractMatrix, theta::Real, opts::RotateOptions, ctx::RunContext)
    out = similar(old_image)
    return prop_rotate!(out, old_image, theta, opts, ctx)
end

@inline function prop_rotate(old_image::AbstractMatrix, theta::Real, opts::RotateOptions)
    out = similar(old_image)
    return prop_rotate!(out, old_image, theta, opts)
end

"""
    prop_rotate(old_image, theta; kwargs...)
    prop_rotate(old_image, theta, ctx; kwargs...)

Return a rotated copy of `old_image`.

# Arguments
- `old_image`: input image
- `theta`: rotation angle in degrees counter-clockwise

# Keywords
- `XC`, `YC`: center of rotation in image pixels
- `XSHIFT`, `YSHIFT`: image shift in pixels
- `MISSING` / `EXTR`: extrapolated fill value
- `METH`: `"linear"` or `"cubic"`
- `CUBIC`: compatibility shortcut for cubic interpolation

# Returns
- An array with the same size and element type as `old_image`.

# Notes
- The default interpolation is linear.
- The cubic option follows the accepted executable-baseline cubic-convolution
  path, not a generic spline interpolation mode.
"""
function prop_rotate(
    old_image::AbstractMatrix,
    theta::Real;
    kwargs...,
)
    opts = isempty(kwargs) ? RotateOptions(old_image) : RotateOptions(old_image, kwargs)
    return prop_rotate(old_image, theta, opts)
end

function prop_rotate(old_image::AbstractMatrix, theta::Real, ctx::RunContext; kwargs...)
    opts = isempty(kwargs) ? RotateOptions(old_image) : RotateOptions(old_image, kwargs)
    return prop_rotate(old_image, theta, opts, ctx)
end
