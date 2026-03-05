@inline function _prop_end_extract_origin(ny::Int, nx::Int, extract::Int)
    extract > 0 || throw(ArgumentError("EXTRACT must be positive"))
    extract <= ny && extract <= nx || throw(ArgumentError("EXTRACT exceeds array size"))
    cy = ny ÷ 2
    cx = nx ÷ 2
    return cy - (extract ÷ 2) + 1, cx - (extract ÷ 2) + 1
end

@inline function _prop_end_layout(ny::Int, nx::Int, extract::Union{Nothing,Int})
    if extract === nothing
        return ny, nx, 1, 1
    end
    r0, c0 = _prop_end_extract_origin(ny, nx, extract)
    return extract, extract, r0, c0
end

@inline function _copy_shifted_complex!(
    out::StridedMatrix{Tout},
    field::StridedMatrix{Tin},
    r0::Int,
    c0::Int,
) where {Tout<:Complex,Tin<:Complex}
    oy, ox = size(out)
    if ka_end_enabled(typeof(out), oy, ox)
        return ka_copy_shifted_complex!(out, field, r0, c0)
    end
    return _copy_shifted_complex_loop!(out, field, r0, c0)
end

@inline function _copy_shifted_complex_loop!(
    out::StridedMatrix{Tout},
    field::StridedMatrix{Tin},
    r0::Int,
    c0::Int,
) where {Tout<:Complex,Tin<:Complex}
    ny, nx = size(field)
    oy, ox = size(out)
    sy = ny ÷ 2
    sx = nx ÷ 2

    @inbounds for j in 1:ox
        js = mod1(c0 + j - 1 + sx, nx)
        for i in 1:oy
            is = mod1(r0 + i - 1 + sy, ny)
            out[i, j] = field[is, js]
        end
    end
    return out
end

@inline function _copy_shifted_intensity!(
    out::StridedMatrix{Tout},
    field::StridedMatrix{Tin},
    r0::Int,
    c0::Int,
) where {Tout<:Number,Tin<:Complex}
    oy, ox = size(out)
    if ka_end_enabled(typeof(out), oy, ox)
        return ka_copy_shifted_intensity!(out, field, r0, c0)
    end
    return _copy_shifted_intensity_loop!(out, field, r0, c0)
end

@inline function _copy_shifted_intensity_loop!(
    out::StridedMatrix{Tout},
    field::StridedMatrix{Tin},
    r0::Int,
    c0::Int,
) where {Tout<:Number,Tin<:Complex}
    ny, nx = size(field)
    oy, ox = size(out)
    sy = ny ÷ 2
    sx = nx ÷ 2

    @inbounds for j in 1:ox
        js = mod1(c0 + j - 1 + sx, nx)
        for i in 1:oy
            is = mod1(r0 + i - 1 + sy, ny)
            out[i, j] = abs2(field[is, js])
        end
    end
    return out
end

@inline function _copy_shifted_complex_generic!(
    out::AbstractMatrix,
    field::AbstractMatrix{<:Complex},
    r0::Int,
    c0::Int,
)
    shifted = prop_shift_center(field)
    @views out .= shifted[r0:(r0 + size(out, 1) - 1), c0:(c0 + size(out, 2) - 1)]
    return out
end

@inline function _copy_shifted_intensity_generic!(
    out::AbstractMatrix,
    field::AbstractMatrix{<:Complex},
    r0::Int,
    c0::Int,
)
    shifted = prop_shift_center(abs2.(field))
    @views out .= shifted[r0:(r0 + size(out, 1) - 1), c0:(c0 + size(out, 2) - 1)]
    return out
end

"""Finalize propagation into a preallocated output buffer."""
function prop_end!(
    out::AbstractMatrix,
    wf::WaveFront;
    noabs::Bool=false,
    extract::Union{Nothing,Int}=nothing,
)
    ny, nx = size(wf.field)
    oy, ox, r0, c0 = _prop_end_layout(ny, nx, extract)
    size(out) == (oy, ox) || throw(ArgumentError("output size must be ($(oy), $(ox))"))

    if noabs
        eltype(out) <: Complex || throw(ArgumentError("output eltype must be complex when noabs=true"))
        if out isa StridedMatrix && wf.field isa StridedMatrix
            return _copy_shifted_complex!(out, wf.field, r0, c0)
        end
        return _copy_shifted_complex_generic!(out, wf.field, r0, c0)
    end

    eltype(out) <: Number || throw(ArgumentError("output eltype must be numeric"))
    if out isa StridedMatrix && wf.field isa StridedMatrix
        return _copy_shifted_intensity!(out, wf.field, r0, c0)
    end
    return _copy_shifted_intensity_generic!(out, wf.field, r0, c0)
end

"""Finalize propagation into caller-provided buffer and return output plus sampling."""
function prop_end(
    wf::WaveFront,
    out::AbstractMatrix;
    noabs::Bool=false,
    extract::Union{Nothing,Int}=nothing,
)
    prop_end!(out, wf; noabs=noabs, extract=extract)
    return out, wf.sampling_m
end

"""Finalize propagation and return either intensity or complex field plus sampling."""
function prop_end(wf::WaveFront; noabs::Bool=false, extract::Union{Nothing,Int}=nothing)
    ny, nx = size(wf.field)
    oy, ox, _, _ = _prop_end_layout(ny, nx, extract)
    out = if noabs
        similar(wf.field, oy, ox)
    else
        RT = typeof(abs2(zero(eltype(wf.field))))
        similar(wf.field, RT, oy, ox)
    end
    prop_end!(out, wf; noabs=noabs, extract=extract)
    return out, wf.sampling_m
end
