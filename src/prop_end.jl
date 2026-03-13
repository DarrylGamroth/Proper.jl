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

@inline function _prop_end_full_quadrants(ny::Int, nx::Int)
    sy = ny ÷ 2
    sx = nx ÷ 2
    return sy, sx, ny - sy, nx - sx
end

abstract type EndFullCopyExecStyle end
struct EndQuadrantCopyStyle <: EndFullCopyExecStyle end
struct EndKernelCopyStyle <: EndFullCopyExecStyle end
struct EndGenericCopyStyle <: EndFullCopyExecStyle end

@inline end_full_copy_exec_style(::Val{false}, ::BackendStyle, ::ArrayLayoutStyle, ::ArrayLayoutStyle) = EndGenericCopyStyle()
@inline end_full_copy_exec_style(::Val{true}, ::CPUBackend, ::StridedLayout, ::StridedLayout) = EndQuadrantCopyStyle()
@inline end_full_copy_exec_style(::Val{true}, ::CUDABackend, ::ArrayLayoutStyle, ::ArrayLayoutStyle) = EndKernelCopyStyle()
@inline end_full_copy_exec_style(::Val{true}, ::BackendStyle, ::ArrayLayoutStyle, ::ArrayLayoutStyle) = EndGenericCopyStyle()

@inline function _copy_shifted_complex_full!(
    out::AbstractMatrix{Tout},
    field::AbstractMatrix{Tin},
) where {Tout<:Complex,Tin<:Complex}
    ny, nx = size(field)
    size(out) == (ny, nx) || throw(ArgumentError("full-copy path requires matching output size"))
    sy, sx, top_h, left_w = _prop_end_full_quadrants(ny, nx)

    @views begin
        copyto!(out[1:top_h, 1:left_w], field[(sy + 1):ny, (sx + 1):nx])
        sx == 0 || copyto!(out[1:top_h, (left_w + 1):nx], field[(sy + 1):ny, 1:sx])
        sy == 0 || copyto!(out[(top_h + 1):ny, 1:left_w], field[1:sy, (sx + 1):nx])
        (sy == 0 || sx == 0) || copyto!(out[(top_h + 1):ny, (left_w + 1):nx], field[1:sy, 1:sx])
    end

    return out
end

@inline function _copy_shifted_intensity_full!(
    out::AbstractMatrix{Tout},
    field::AbstractMatrix{Tin},
) where {Tout<:Number,Tin<:Complex}
    ny, nx = size(field)
    size(out) == (ny, nx) || throw(ArgumentError("full-copy path requires matching output size"))
    sy, sx, top_h, left_w = _prop_end_full_quadrants(ny, nx)

    @views begin
        broadcast!(abs2, out[1:top_h, 1:left_w], field[(sy + 1):ny, (sx + 1):nx])
        sx == 0 || broadcast!(abs2, out[1:top_h, (left_w + 1):nx], field[(sy + 1):ny, 1:sx])
        sy == 0 || broadcast!(abs2, out[(top_h + 1):ny, 1:left_w], field[1:sy, (sx + 1):nx])
        (sy == 0 || sx == 0) || broadcast!(abs2, out[(top_h + 1):ny, (left_w + 1):nx], field[1:sy, 1:sx])
    end

    return out
end

@inline function end_full_copy_exec_style(::Type{A}, ::Type{B}) where {A<:AbstractArray,B<:AbstractArray}
    return end_full_copy_exec_style(
        Val(same_backend_style(A, B)),
        backend_style(A),
        array_layout_style(A),
        array_layout_style(B),
    )
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
    @views copyto!(out, shifted[r0:(r0 + size(out, 1) - 1), c0:(c0 + size(out, 2) - 1)])
    return out
end

@inline function _copy_shifted_intensity_generic!(
    out::AbstractMatrix,
    field::AbstractMatrix{<:Complex},
    r0::Int,
    c0::Int,
)
    shifted = prop_shift_center(abs2.(field))
    @views copyto!(out, shifted[r0:(r0 + size(out, 1) - 1), c0:(c0 + size(out, 2) - 1)])
    return out
end

@inline function _prop_end_copy_noabs!(
    out::StridedMatrix{<:Complex},
    field::StridedMatrix{<:Complex},
    r0::Int,
    c0::Int,
)
    return _copy_shifted_complex!(out, field, r0, c0)
end

@inline _prop_end_full_copy_noabs!(::EndQuadrantCopyStyle, out::AbstractMatrix{<:Complex}, field::AbstractMatrix{<:Complex}) =
    _copy_shifted_complex_full!(out, field)

@inline _prop_end_full_copy_noabs!(::EndKernelCopyStyle, out::AbstractMatrix{<:Complex}, field::AbstractMatrix{<:Complex}) =
    _copy_shifted_complex!(out, field, 1, 1)

@inline _prop_end_full_copy_noabs!(::EndGenericCopyStyle, out::AbstractMatrix{<:Complex}, field::AbstractMatrix{<:Complex}) =
    _copy_shifted_complex_generic!(out, field, 1, 1)

@inline function _prop_end_copy_noabs!(
    out::AbstractMatrix,
    field::AbstractMatrix{<:Complex},
    r0::Int,
    c0::Int,
)
    oy, ox = size(out)
    if ka_end_enabled(typeof(out), oy, ox) && same_backend_style(typeof(out), typeof(field))
        return ka_copy_shifted_complex!(out, field, r0, c0)
    end
    return _copy_shifted_complex_generic!(out, field, r0, c0)
end

@inline function _prop_end_copy_intensity!(
    out::StridedMatrix,
    field::StridedMatrix{<:Complex},
    r0::Int,
    c0::Int,
)
    return _copy_shifted_intensity!(out, field, r0, c0)
end

@inline _prop_end_full_copy_intensity!(::EndQuadrantCopyStyle, out::AbstractMatrix, field::AbstractMatrix{<:Complex}) =
    _copy_shifted_intensity_full!(out, field)

@inline _prop_end_full_copy_intensity!(::EndKernelCopyStyle, out::AbstractMatrix, field::AbstractMatrix{<:Complex}) =
    _copy_shifted_intensity!(out, field, 1, 1)

@inline _prop_end_full_copy_intensity!(::EndGenericCopyStyle, out::AbstractMatrix, field::AbstractMatrix{<:Complex}) =
    _copy_shifted_intensity_generic!(out, field, 1, 1)

@inline function _prop_end_copy_intensity!(
    out::AbstractMatrix,
    field::AbstractMatrix{<:Complex},
    r0::Int,
    c0::Int,
)
    oy, ox = size(out)
    if ka_end_enabled(typeof(out), oy, ox) && same_backend_style(typeof(out), typeof(field))
        return ka_copy_shifted_intensity!(out, field, r0, c0)
    end
    return _copy_shifted_intensity_generic!(out, field, r0, c0)
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
    fullcopy_sty = end_full_copy_exec_style(typeof(out), typeof(wf.field))
    fullcopy = extract === nothing

    if noabs
        eltype(out) <: Complex || throw(ArgumentError("output eltype must be complex when noabs=true"))
        fullcopy && return _prop_end_full_copy_noabs!(fullcopy_sty, out, wf.field)
        return _prop_end_copy_noabs!(out, wf.field, r0, c0)
    end

    eltype(out) <: Number || throw(ArgumentError("output eltype must be numeric"))
    fullcopy && return _prop_end_full_copy_intensity!(fullcopy_sty, out, wf.field)
    return _prop_end_copy_intensity!(out, wf.field, r0, c0)
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
