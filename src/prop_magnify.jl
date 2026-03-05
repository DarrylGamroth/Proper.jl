"""Magnify image with cubic-convolution interpolation and optional flux scaling."""
function prop_magnify(image_in::AbstractMatrix, mag0::Real, size_out0::Integer=0; kwargs...)
    mag = float(mag0)
    ny, nx = size(image_in)
    out_n = size_out0 > 0 ? Int(size_out0) : round(Int, ny * mag)
    out_n > 0 || throw(ArgumentError("size_out must be positive"))

    Tin = typeof(real(zero(eltype(image_in))))
    T = float(promote_type(typeof(mag), Tin))
    cx_in = T(nx ÷ 2) + one(T)
    cy_in = T(ny ÷ 2) + one(T)
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

    out = prop_cubic_conv(image_in, xcoords, ycoords; grid=true)

    if switch_set(:CONSERVE; kwargs...)
        if eltype(image_in) <: Complex
            out ./= mag
        else
            out ./= mag^2
        end
    elseif switch_set(:AMP_CONSERVE; kwargs...)
        out ./= mag
    end

    return out
end
