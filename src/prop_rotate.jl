"""Rotate image by theta degrees around center using bilinear interpolation."""
function prop_rotate(old_image::AbstractMatrix, theta::Real; kwargs...)
    ny, nx = size(old_image)
    out = similar(old_image)

    ang = deg2rad(float(theta))
    c = cos(ang)
    s = sin(ang)
    cy = (ny + 1) / 2
    cx = (nx + 1) / 2

    @inbounds for j in 1:nx
        x = j - cx
        for i in 1:ny
            y = i - cy
            xr = c * x + s * y + cx
            yr = -s * x + c * y + cy
            out[i, j] = bilinear_sample(old_image, yr, xr)
        end
    end

    return out
end
