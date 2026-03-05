"""Cubic convolution interpolation over scalar, axis-grid, or coordinate-grid requests."""
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

function prop_cubic_conv_grid!(out::StridedMatrix, sty::InterpStyle, a::StridedMatrix, xval::AbstractVector, yval::AbstractVector)
    size(out) == (length(yval), length(xval)) || throw(ArgumentError("output size mismatch for grid interpolation"))
    @inbounds for j in eachindex(xval)
        x = xval[j]
        for i in eachindex(yval)
            out[i, j] = _cubic_sample(sty, a, yval[i], x)
        end
    end
    return out
end

function prop_cubic_conv_grid!(out::AbstractMatrix, sty::InterpStyle, a::AbstractMatrix, xval::AbstractVector, yval::AbstractVector)
    size(out) == (length(yval), length(xval)) || throw(ArgumentError("output size mismatch for grid interpolation"))
    if out isa StridedMatrix && a isa StridedMatrix
        return prop_cubic_conv_grid!(out, sty, a, xval, yval)
    end
    host_out = Matrix{eltype(out)}(undef, size(out)...)
    prop_cubic_conv_grid!(host_out, sty, Matrix(a), xval, yval)
    copyto!(out, host_out)
    return out
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

function prop_cubic_conv(sty::InterpStyle, a::AbstractMatrix, y::Real, x::Real)
    return a isa StridedMatrix ? _cubic_sample(sty, a, y, x) : _cubic_sample(sty, Matrix(a), y, x)
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
    if a isa StridedMatrix
        return _prop_cubic_conv(sty, PointwiseTopology(), a, xval, yval)
    end
    host_out = _prop_cubic_conv(sty, PointwiseTopology(), Matrix(a), xval, yval)
    out = similar(a, eltype(host_out), size(host_out)...)
    copyto!(out, host_out)
    return out
end

function prop_cubic_conv(a::AbstractMatrix, xval::AbstractVector, yval::AbstractVector; threaded::Bool=true, grid::Bool=true)
    return prop_cubic_conv(_default_interp_style(a), a, xval, yval; threaded=threaded, grid=grid)
end

function prop_cubic_conv(ctx::RunContext, a::AbstractMatrix, xval::AbstractVector, yval::AbstractVector; threaded::Bool=true, grid::Bool=true)
    return prop_cubic_conv(interp_style(ctx), a, xval, yval; threaded=threaded, grid=grid)
end

function prop_cubic_conv(sty::InterpStyle, a::AbstractMatrix, xgrid::AbstractMatrix, ygrid::AbstractMatrix; threaded::Bool=true, grid::Bool=false)
    if a isa StridedMatrix
        return _prop_cubic_conv(sty, PointwiseTopology(), a, xgrid, ygrid)
    end
    host_out = _prop_cubic_conv(sty, PointwiseTopology(), Matrix(a), xgrid, ygrid)
    out = similar(a, eltype(host_out), size(host_out)...)
    copyto!(out, host_out)
    return out
end

function prop_cubic_conv(a::AbstractMatrix, xgrid::AbstractMatrix, ygrid::AbstractMatrix; threaded::Bool=true, grid::Bool=false)
    return prop_cubic_conv(_default_interp_style(a), a, xgrid, ygrid; threaded=threaded, grid=grid)
end

function prop_cubic_conv(ctx::RunContext, a::AbstractMatrix, xgrid::AbstractMatrix, ygrid::AbstractMatrix; threaded::Bool=true, grid::Bool=false)
    return prop_cubic_conv(interp_style(ctx), a, xgrid, ygrid; threaded=threaded, grid=grid)
end
