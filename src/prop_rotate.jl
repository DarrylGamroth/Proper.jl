@inline function _rotate_method(kwargs)
    if haskey(kwargs, :METH)
        return Symbol(lowercase(String(kwargs[:METH])))
    elseif haskey(kwargs, :meth)
        return Symbol(lowercase(String(kwargs[:meth])))
    elseif haskey(kwargs, :CUBIC) || haskey(kwargs, :cubic)
        return :cubic
    end
    return :cubic
end

"""Rotate image by theta degrees around center (cubic by default)."""
function prop_rotate(old_image::AbstractMatrix, theta::Real; kwargs...)
    ny, nx = size(old_image)
    out = similar(old_image)

    meth = _rotate_method(kwargs)
    ang = deg2rad(-float(theta))
    c = cos(ang)
    s = sin(ang)
    cx = haskey(kwargs, :XC) ? float(kwargs[:XC]) : haskey(kwargs, :xc) ? float(kwargs[:xc]) : (nx ÷ 2 + 1)
    cy = haskey(kwargs, :YC) ? float(kwargs[:YC]) : haskey(kwargs, :yc) ? float(kwargs[:yc]) : (ny ÷ 2 + 1)
    sx = haskey(kwargs, :XSHIFT) ? float(kwargs[:XSHIFT]) : haskey(kwargs, :xshift) ? float(kwargs[:xshift]) : 0.0
    sy = haskey(kwargs, :YSHIFT) ? float(kwargs[:YSHIFT]) : haskey(kwargs, :yshift) ? float(kwargs[:yshift]) : 0.0

    @inbounds for j in 1:nx
        x = j - cx - sx
        for i in 1:ny
            y = i - cy - sy
            xr = c * x - s * y + cx
            yr = s * x + c * y + cy
            if meth === :linear
                out[i, j] = bilinear_sample(old_image, yr, xr)
            else
                out[i, j] = prop_cubic_conv(old_image, yr, xr)
            end
        end
    end

    return out
end
