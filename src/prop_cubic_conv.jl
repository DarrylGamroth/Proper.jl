"""Cubic convolution interpolation over scalar, axis-grid, or coordinate-grid requests."""
function prop_cubic_conv(a::AbstractMatrix, y::Real, x::Real)
    return libcconv(a, y, x)
end

function prop_cubic_conv(a::AbstractMatrix, xval::AbstractVector, yval::AbstractVector; threaded::Bool=true, grid::Bool=true)
    if !grid
        length(xval) == length(yval) || throw(ArgumentError("xval and yval lengths must match when grid=false"))
        out = similar(a, length(yval))
        @inbounds for i in eachindex(xval)
            out[i] = libcconv(a, yval[i], xval[i])
        end
        return out
    end

    out = similar(a, length(yval), length(xval))
    @inbounds for j in eachindex(xval)
        x = xval[j]
        for i in eachindex(yval)
            out[i, j] = libcconv(a, yval[i], x)
        end
    end
    return out
end

function prop_cubic_conv(a::AbstractMatrix, xgrid::AbstractMatrix, ygrid::AbstractMatrix; threaded::Bool=true, grid::Bool=false)
    size(xgrid) == size(ygrid) || throw(ArgumentError("xgrid and ygrid sizes must match"))
    out = similar(a, size(xgrid)...)
    @inbounds for j in axes(xgrid, 2)
        for i in axes(xgrid, 1)
            out[i, j] = libcconv(a, ygrid[i, j], xgrid[i, j])
        end
    end
    return out
end
