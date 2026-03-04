@inline function _clamp_index(i::Int, n::Int)
    return ifelse(i < 1, 1, ifelse(i > n, n, i))
end

@inline function bilinear_sample(a::AbstractMatrix{T}, y::Real, x::Real) where {T}
    ny, nx = size(a)
    x0 = floor(Int, x)
    y0 = floor(Int, y)
    x1 = x0 + 1
    y1 = y0 + 1

    x0c = _clamp_index(x0, nx)
    x1c = _clamp_index(x1, nx)
    y0c = _clamp_index(y0, ny)
    y1c = _clamp_index(y1, ny)

    tx = x - x0
    ty = y - y0

    v00 = a[y0c, x0c]
    v10 = a[y0c, x1c]
    v01 = a[y1c, x0c]
    v11 = a[y1c, x1c]

    v0 = (1 - tx) * v00 + tx * v10
    v1 = (1 - tx) * v01 + tx * v11
    return (1 - ty) * v0 + ty * v1
end
