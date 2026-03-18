abstract type CubicConvTopology end
struct GridTopology <: CubicConvTopology end
struct PointwiseTopology <: CubicConvTopology end

@inline function _cubic_sample(::GenericInterpStyle, a::AbstractMatrix, y::Real, x::Real)
    return libcconv(a, y, x)
end

@inline function _cubic_sample(::CubicInterpStyle, a::AbstractMatrix, y::Real, x::Real)
    return libcconv(a, y, x)
end

@inline _default_interp_style(a::AbstractMatrix) = interp_style(typeof(a))

function _prop_cubic_conv(sty::InterpStyle, ::PointwiseTopology, a::StridedMatrix, xval::AbstractVector, yval::AbstractVector)
    length(xval) == length(yval) || throw(ArgumentError("xval and yval lengths must match when grid=false"))
    out = similar(a, length(yval))
    @inbounds for i in eachindex(xval)
        out[i] = _cubic_sample(sty, a, yval[i], xval[i])
    end
    return out
end

function _prop_cubic_conv_grid_loop!(
    out::StridedMatrix,
    sty::InterpStyle,
    a::StridedMatrix,
    xval::AbstractVector,
    yval::AbstractVector,
)
    size(out) == (length(yval), length(xval)) || throw(ArgumentError("output size mismatch for grid interpolation"))
    @inbounds for j in eachindex(xval)
        x = xval[j]
        for i in eachindex(yval)
            out[i, j] = _cubic_sample(sty, a, yval[i], x)
        end
    end
    return out
end

abstract type CubicGridExecStyle end
struct CubicGridLoopExecStyle <: CubicGridExecStyle end
struct CubicGridKAExecStyle <: CubicGridExecStyle end
struct CubicGridHostExecStyle <: CubicGridExecStyle end
struct CubicGridUnsupportedExecStyle <: CubicGridExecStyle end

@inline cubic_grid_exec_style(
    ::StridedLayout,
    ::StridedLayout,
    ::CPUBackend,
    ::CPUBackend,
    ::InterpStyle,
    ::Val{false},
) = CubicGridLoopExecStyle()

@inline cubic_grid_exec_style(
    ::ArrayLayoutStyle,
    ::ArrayLayoutStyle,
    ::B,
    ::B,
    ::CubicInterpStyle,
    ::Val{true},
) where {B<:BackendStyle} = CubicGridKAExecStyle()

@inline cubic_grid_exec_style(
    ::ArrayLayoutStyle,
    ::ArrayLayoutStyle,
    ::CPUBackend,
    ::CPUBackend,
    ::InterpStyle,
    ::Val,
) = CubicGridHostExecStyle()

@inline cubic_grid_exec_style(
    ::ArrayLayoutStyle,
    ::ArrayLayoutStyle,
    ::BackendStyle,
    ::BackendStyle,
    ::InterpStyle,
    ::Val,
) = CubicGridUnsupportedExecStyle()

@inline function _prop_cubic_conv_grid!(
    ::CubicGridLoopExecStyle,
    out::StridedMatrix,
    sty::InterpStyle,
    a::StridedMatrix,
    xval::AbstractVector,
    yval::AbstractVector,
)
    return _prop_cubic_conv_grid_loop!(out, sty, a, xval, yval)
end

@inline function _prop_cubic_conv_grid!(
    ::CubicGridKAExecStyle,
    out::AbstractMatrix,
    ::CubicInterpStyle,
    a::AbstractMatrix,
    xval::AbstractVector,
    yval::AbstractVector,
)
    return ka_cubic_conv_grid!(out, a, backend_adapt(out, xval), backend_adapt(out, yval))
end

@inline function _prop_cubic_conv_grid!(
    ::CubicGridHostExecStyle,
    out::AbstractMatrix,
    sty::InterpStyle,
    a::AbstractMatrix,
    xval::AbstractVector,
    yval::AbstractVector,
)
    host_out = Matrix{eltype(out)}(undef, size(out)...)
    prop_cubic_conv_grid!(host_out, sty, Matrix(a), xval, yval)
    copyto!(out, host_out)
    return out
end

@inline function _prop_cubic_conv_grid!(
    ::CubicGridUnsupportedExecStyle,
    out::AbstractMatrix,
    sty::InterpStyle,
    a::AbstractMatrix,
    xval::AbstractVector,
    yval::AbstractVector,
)
    _ = sty
    _ = xval
    _ = yval
    throw(ArgumentError("prop_cubic_conv_grid! has no native implementation for $(typeof(a)) -> $(typeof(out)); use matching CPU arrays or a backend-native path"))
end

"""
    prop_cubic_conv_grid!(out, sty, a, xval, yval)
    prop_cubic_conv_grid!(out, a, xval, yval)
    prop_cubic_conv_grid!(out, ctx, a, xval, yval)

Evaluate cubic-convolution interpolation over the tensor-product grid defined
by `xval` and `yval` and write the result into `out`.

# Arguments
- `out`: destination array of size `(length(yval), length(xval))`
- `a`: source image
- `xval`, `yval`: grid axes in the cubic-convolution coordinate system
"""
function prop_cubic_conv_grid!(out::AbstractMatrix, sty::InterpStyle, a::AbstractMatrix, xval::AbstractVector, yval::AbstractVector)
    size(out) == (length(yval), length(xval)) || throw(ArgumentError("output size mismatch for grid interpolation"))
    oy, ox = size(out)
    sty_exec = cubic_grid_exec_style(
        array_layout_style(typeof(out)),
        array_layout_style(typeof(a)),
        backend_style(typeof(out)),
        backend_style(typeof(a)),
        sty,
        Val(ka_cubic_grid_enabled(typeof(out), oy, ox)),
    )
    return _prop_cubic_conv_grid!(sty_exec, out, sty, a, xval, yval)
end

@inline function prop_cubic_conv_grid!(out::AbstractMatrix, a::AbstractMatrix, xval::AbstractVector, yval::AbstractVector)
    return prop_cubic_conv_grid!(out, _default_interp_style(a), a, xval, yval)
end

@inline function prop_cubic_conv_grid!(out::AbstractMatrix, ctx::RunContext, a::AbstractMatrix, xval::AbstractVector, yval::AbstractVector)
    return prop_cubic_conv_grid!(out, interp_style(ctx), a, xval, yval)
end

function _prop_cubic_conv(sty::InterpStyle, ::PointwiseTopology, a::StridedMatrix, xgrid::AbstractMatrix, ygrid::AbstractMatrix)
    size(xgrid) == size(ygrid) || throw(ArgumentError("xgrid and ygrid sizes must match"))
    out = similar(a, size(xgrid)...)
    @inbounds for j in axes(xgrid, 2)
        for i in axes(xgrid, 1)
            out[i, j] = _cubic_sample(sty, a, ygrid[i, j], xgrid[i, j])
        end
    end
    return out
end

"""
    prop_cubic_conv(args...)

Cubic-convolution interpolation over scalar, axis-grid, or coordinate-grid
requests.

# Notes
- This is the public wrapper around the upstream PROPER cubic-convolution
  kernel.
- Scalar calls sample one point, vector calls support pointwise or grid
  interpolation, and `prop_cubic_conv_grid!` writes directly into a caller
  buffer.
- GPU/backend execution requires a native implementation for the requested
  topology. Unsupported backend combinations throw instead of silently falling
  back to host materialization.
"""
function prop_cubic_conv(sty::InterpStyle, a::AbstractMatrix, y::Real, x::Real)
    return _prop_cubic_conv_point(array_layout_style(typeof(a)), backend_style(typeof(a)), sty, a, y, x)
end

@inline function _prop_cubic_conv_point(
    ::StridedLayout,
    ::CPUBackend,
    sty::InterpStyle,
    a::StridedMatrix,
    y::Real,
    x::Real,
)
    return _cubic_sample(sty, a, y, x)
end

@inline function _prop_cubic_conv_point(
    ::ArrayLayoutStyle,
    ::CPUBackend,
    sty::InterpStyle,
    a::AbstractMatrix,
    y::Real,
    x::Real,
)
    return _cubic_sample(sty, Matrix(a), y, x)
end

@inline function _prop_cubic_conv_point(
    ::ArrayLayoutStyle,
    ::BackendStyle,
    sty::InterpStyle,
    a::AbstractMatrix,
    y::Real,
    x::Real,
)
    _ = sty
    _ = y
    _ = x
    throw(ArgumentError("prop_cubic_conv has no native point-sampling implementation for $(typeof(a)); use CPU arrays or a backend-native path"))
end

function prop_cubic_conv(a::AbstractMatrix, y::Real, x::Real)
    return prop_cubic_conv(_default_interp_style(a), a, y, x)
end

function prop_cubic_conv(ctx::RunContext, a::AbstractMatrix, y::Real, x::Real)
    return prop_cubic_conv(interp_style(ctx), a, y, x)
end

function prop_cubic_conv(sty::InterpStyle, a::AbstractMatrix, xval::AbstractVector, yval::AbstractVector; threaded::Bool=true, grid::Bool=true)
    if grid
        out = similar(a, length(yval), length(xval))
        return prop_cubic_conv_grid!(out, sty, a, xval, yval)
    end
    return _prop_cubic_conv_pointwise(array_layout_style(typeof(a)), backend_style(typeof(a)), sty, a, xval, yval)
end

@inline function _prop_cubic_conv_pointwise(
    ::StridedLayout,
    ::CPUBackend,
    sty::InterpStyle,
    a::StridedMatrix,
    xval::AbstractVector,
    yval::AbstractVector,
)
    return _prop_cubic_conv(sty, PointwiseTopology(), a, xval, yval)
end

@inline function _prop_cubic_conv_pointwise(
    ::ArrayLayoutStyle,
    ::CPUBackend,
    sty::InterpStyle,
    a::AbstractMatrix,
    xval::AbstractVector,
    yval::AbstractVector,
)
    host_out = _prop_cubic_conv(sty, PointwiseTopology(), Matrix(a), xval, yval)
    out = similar(a, eltype(host_out), size(host_out)...)
    copyto!(out, host_out)
    return out
end

@inline function _prop_cubic_conv_pointwise(
    ::ArrayLayoutStyle,
    ::BackendStyle,
    sty::InterpStyle,
    a::AbstractMatrix,
    xval::AbstractVector,
    yval::AbstractVector,
)
    _ = sty
    _ = xval
    _ = yval
    throw(ArgumentError("prop_cubic_conv pointwise mode has no native implementation for $(typeof(a)); use CPU arrays or a backend-native path"))
end

function prop_cubic_conv(a::AbstractMatrix, xval::AbstractVector, yval::AbstractVector; threaded::Bool=true, grid::Bool=true)
    return prop_cubic_conv(_default_interp_style(a), a, xval, yval; threaded=threaded, grid=grid)
end

function prop_cubic_conv(ctx::RunContext, a::AbstractMatrix, xval::AbstractVector, yval::AbstractVector; threaded::Bool=true, grid::Bool=true)
    return prop_cubic_conv(interp_style(ctx), a, xval, yval; threaded=threaded, grid=grid)
end

function prop_cubic_conv(sty::InterpStyle, a::AbstractMatrix, xgrid::AbstractMatrix, ygrid::AbstractMatrix; threaded::Bool=true, grid::Bool=false)
    return _prop_cubic_conv_coordinate_grid(array_layout_style(typeof(a)), backend_style(typeof(a)), sty, a, xgrid, ygrid)
end

@inline function _prop_cubic_conv_coordinate_grid(
    ::StridedLayout,
    ::CPUBackend,
    sty::InterpStyle,
    a::StridedMatrix,
    xgrid::AbstractMatrix,
    ygrid::AbstractMatrix,
)
    return _prop_cubic_conv(sty, PointwiseTopology(), a, xgrid, ygrid)
end

@inline function _prop_cubic_conv_coordinate_grid(
    ::ArrayLayoutStyle,
    ::CPUBackend,
    sty::InterpStyle,
    a::AbstractMatrix,
    xgrid::AbstractMatrix,
    ygrid::AbstractMatrix,
)
    host_out = _prop_cubic_conv(sty, PointwiseTopology(), Matrix(a), xgrid, ygrid)
    out = similar(a, eltype(host_out), size(host_out)...)
    copyto!(out, host_out)
    return out
end

@inline function _prop_cubic_conv_coordinate_grid(
    ::ArrayLayoutStyle,
    ::BackendStyle,
    sty::InterpStyle,
    a::AbstractMatrix,
    xgrid::AbstractMatrix,
    ygrid::AbstractMatrix,
)
    _ = sty
    _ = xgrid
    _ = ygrid
    throw(ArgumentError("prop_cubic_conv grid=false coordinate mode has no native implementation for $(typeof(a)); use CPU arrays or a backend-native path"))
end

function prop_cubic_conv(a::AbstractMatrix, xgrid::AbstractMatrix, ygrid::AbstractMatrix; threaded::Bool=true, grid::Bool=false)
    return prop_cubic_conv(_default_interp_style(a), a, xgrid, ygrid; threaded=threaded, grid=grid)
end

function prop_cubic_conv(ctx::RunContext, a::AbstractMatrix, xgrid::AbstractMatrix, ygrid::AbstractMatrix; threaded::Bool=true, grid::Bool=false)
    return prop_cubic_conv(interp_style(ctx), a, xgrid, ygrid; threaded=threaded, grid=grid)
end
