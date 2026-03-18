@inline _center_shift_amount(n::Integer, inverse::Bool) = inverse ? cld(n, 2) : fld(n, 2)

@inline function _shift_center_exec_style(::Type{A}, ::Type{B}, ny::Integer, nx::Integer) where {A<:AbstractArray,B<:AbstractArray}
    return ka_shift_enabled(A, ny, nx) && same_backend_style(A, B) ? ShiftKAStyle() : ShiftLoopStyle()
end

"""Circularly shift array from origin to center, or back with `inverse=true`."""
function prop_shift_center(a::AbstractMatrix; inverse::Bool=false)
    out = similar(a, size(a))
    return prop_shift_center!(out, a; inverse=inverse)
end

function prop_shift_center!(out::AbstractMatrix{T}, a::AbstractMatrix; inverse::Bool=false) where {T}
    size(out) == size(a) || throw(ArgumentError("output size must match input size"))
    sy = _center_shift_amount(size(a, 1), inverse)
    sx = _center_shift_amount(size(a, 2), inverse)
    return _prop_shift_center!(_shift_center_exec_style(typeof(out), typeof(a), size(a, 1), size(a, 2)), out, a, sy, sx)
end

"""
    prop_shift_center!(out, a; inverse=false)

Circularly shift `a` between origin ordering and centered ordering, writing the
result to `out`.

# Arguments
- `out`: destination array with the same size as `a`
- `a`: input array

# Keywords
- `inverse`: when `true`, shift from centered ordering back to origin ordering
"""
function prop_shift_center!(out::StridedMatrix{T}, a::StridedMatrix; inverse::Bool=false) where {T}
    size(out) == size(a) || throw(ArgumentError("output size must match input size"))
    sy = _center_shift_amount(size(a, 1), inverse)
    sx = _center_shift_amount(size(a, 2), inverse)
    return shift_copy!(out, a, sy, sx)
end

@inline function _prop_shift_center!(::ShiftLoopStyle, out::AbstractMatrix{T}, a::AbstractMatrix, sy::Integer, sx::Integer) where {T}
    @inbounds for j in axes(a, 2)
        js = mod1(j + sx, size(a, 2))
        for i in axes(a, 1)
            is = mod1(i + sy, size(a, 1))
            out[is, js] = a[i, j]
        end
    end
    return out
end

@inline _prop_shift_center!(::ShiftKAStyle, out::AbstractMatrix, a::AbstractMatrix, sy::Integer, sx::Integer) =
    ka_shift_copy!(out, a, sy, sx)

@inline function shift_center_for_wavefront!(wf::WaveFront, a::AbstractMatrix; inverse::Bool=false)
    same_backend_style(typeof(wf.field), typeof(a)) || return backend_adapt(wf.field, prop_shift_center(a; inverse=inverse))
    scratch = if eltype(a) <: Real
        ensure_fft_real_scratch!(wf.workspace.fft, size(a, 1), size(a, 2))
    else
        ensure_fft_scratch!(wf.workspace.fft, size(a, 1), size(a, 2))
    end
    prop_shift_center!(scratch, a; inverse=inverse)
    return scratch
end
