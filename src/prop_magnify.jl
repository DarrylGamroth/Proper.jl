"""Bilinear image magnification with optional output size."""
function prop_magnify(image_in::AbstractMatrix, mag0::Real, size_out0::Integer=0; kwargs...)
    mag = float(mag0)
    ny, nx = size(image_in)
    out_n = size_out0 > 0 ? Int(size_out0) : round(Int, ny * mag)
    out = similar(image_in, out_n, out_n)

    cy_in = (ny + 1) / 2
    cx_in = (nx + 1) / 2
    cy_out = (out_n + 1) / 2
    cx_out = (out_n + 1) / 2

    @inbounds for j in 1:out_n
        x_in = (j - cx_out) / mag + cx_in
        for i in 1:out_n
            y_in = (i - cy_out) / mag + cy_in
            out[i, j] = bilinear_sample(image_in, y_in, x_in)
        end
    end
    return out
end
