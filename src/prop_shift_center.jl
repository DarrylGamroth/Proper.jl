@inline _center_shift_amount(n::Integer, inverse::Bool) = inverse ? cld(n, 2) : fld(n, 2)

"""Circularly shift array from origin to center, or back with `inverse=true`."""
function prop_shift_center(a::AbstractMatrix; inverse::Bool=false)
    out = similar(a, size(a))
    return prop_shift_center!(out, a; inverse=inverse)
end

function prop_shift_center!(out::AbstractMatrix{T}, a::AbstractMatrix; inverse::Bool=false) where {T}
    size(out) == size(a) || throw(ArgumentError("output size must match input size"))
    sy = _center_shift_amount(size(a, 1), inverse)
    sx = _center_shift_amount(size(a, 2), inverse)
    @inbounds for j in axes(a, 2)
        js = mod1(j + sx, size(a, 2))
        for i in axes(a, 1)
            is = mod1(i + sy, size(a, 1))
            out[is, js] = a[i, j]
        end
    end
    return out
end

function prop_shift_center!(out::StridedMatrix{T}, a::StridedMatrix; inverse::Bool=false) where {T}
    size(out) == size(a) || throw(ArgumentError("output size must match input size"))
    sy = _center_shift_amount(size(a, 1), inverse)
    sx = _center_shift_amount(size(a, 2), inverse)
    return shift_copy!(out, a, sy, sx)
end
