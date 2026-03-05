"""Cubic convolution interpolation over scalar, axis-grid, or coordinate-grid requests."""
abstract type CubicConvTopology end
struct GridTopology <: CubicConvTopology end
struct PointwiseTopology <: CubicConvTopology end

@inline _interp_style_for(::Type{A}) where {A<:AbstractMatrix} = interp_style(A)

@inline function _cubic_sample(::GenericInterpStyle, a::AbstractMatrix, y::Real, x::Real)
    return libcconv(a, y, x)
end

@inline function _cubic_sample(::CubicInterpStyle, a::AbstractMatrix, y::Real, x::Real)
    return libcconv(a, y, x)
end

function _prop_cubic_conv(::InterpStyle, ::PointwiseTopology, a::AbstractMatrix, xval::AbstractVector, yval::AbstractVector)
    length(xval) == length(yval) || throw(ArgumentError("xval and yval lengths must match when grid=false"))
    out = similar(a, length(yval))
    sty = _interp_style_for(typeof(a))
    @inbounds for i in eachindex(xval)
        out[i] = _cubic_sample(sty, a, yval[i], xval[i])
    end
    return out
end

function _prop_cubic_conv(::InterpStyle, ::GridTopology, a::AbstractMatrix, xval::AbstractVector, yval::AbstractVector)
    out = similar(a, length(yval), length(xval))
    sty = _interp_style_for(typeof(a))
    @inbounds for j in eachindex(xval)
        x = xval[j]
        for i in eachindex(yval)
            out[i, j] = _cubic_sample(sty, a, yval[i], x)
        end
    end
    return out
end

function _prop_cubic_conv(::InterpStyle, ::PointwiseTopology, a::AbstractMatrix, xgrid::AbstractMatrix, ygrid::AbstractMatrix)
    size(xgrid) == size(ygrid) || throw(ArgumentError("xgrid and ygrid sizes must match"))
    out = similar(a, size(xgrid)...)
    sty = _interp_style_for(typeof(a))
    @inbounds for j in axes(xgrid, 2)
        for i in axes(xgrid, 1)
            out[i, j] = _cubic_sample(sty, a, ygrid[i, j], xgrid[i, j])
        end
    end
    return out
end

function prop_cubic_conv(a::AbstractMatrix, y::Real, x::Real)
    return _cubic_sample(_interp_style_for(typeof(a)), a, y, x)
end

function prop_cubic_conv(a::AbstractMatrix, xval::AbstractVector, yval::AbstractVector; threaded::Bool=true, grid::Bool=true)
    topo = grid ? GridTopology() : PointwiseTopology()
    return _prop_cubic_conv(_interp_style_for(typeof(a)), topo, a, xval, yval)
end

function prop_cubic_conv(a::AbstractMatrix, xgrid::AbstractMatrix, ygrid::AbstractMatrix; threaded::Bool=true, grid::Bool=false)
    return _prop_cubic_conv(_interp_style_for(typeof(a)), PointwiseTopology(), a, xgrid, ygrid)
end
