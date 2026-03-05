using KernelAbstractions
import AcceleratedKernels as AK

@kernel function _ka_apply_shifted_mask_kernel!(
    field,
    mask,
    sy::Int,
    sx::Int,
    ny::Int,
    nx::Int,
    invert::Bool,
)
    I = @index(Global, NTuple)
    i = I[1]
    j = I[2]

    if i <= ny && j <= nx
        is = i + sy
        js = j + sx
        is = ifelse(is > ny, is - ny, is)
        js = ifelse(js > nx, js - nx, js)
        m = mask[is, js]
        field[i, j] *= invert ? (one(m) - m) : m
    end
end

@kernel function _ka_copy_shifted_complex_kernel!(
    out,
    field,
    r0::Int,
    c0::Int,
    sy::Int,
    sx::Int,
    ny::Int,
    nx::Int,
    oy::Int,
    ox::Int,
)
    I = @index(Global, NTuple)
    i = I[1]
    j = I[2]

    if i <= oy && j <= ox
        is = r0 + i - 1 + sy
        js = c0 + j - 1 + sx
        is = ifelse(is > ny, is - ny, is)
        js = ifelse(js > nx, js - nx, js)
        out[i, j] = field[is, js]
    end
end

@kernel function _ka_copy_shifted_intensity_kernel!(
    out,
    field,
    r0::Int,
    c0::Int,
    sy::Int,
    sx::Int,
    ny::Int,
    nx::Int,
    oy::Int,
    ox::Int,
)
    I = @index(Global, NTuple)
    i = I[1]
    j = I[2]

    if i <= oy && j <= ox
        is = r0 + i - 1 + sy
        js = c0 + j - 1 + sx
        is = ifelse(is > ny, is - ny, is)
        js = ifelse(js > nx, js - nx, js)
        out[i, j] = abs2(field[is, js])
    end
end

@inline function ka_apply_shifted_mask!(
    field::AbstractMatrix{<:Complex},
    mask::AbstractMatrix{<:Real};
    invert::Bool=false,
)
    ny, nx = size(field)
    backend = AK.get_backend(field)
    _ka_apply_shifted_mask_kernel!(backend, (16, 16))(field, mask, ny ÷ 2, nx ÷ 2, ny, nx, invert; ndrange=(ny, nx))
    AK.synchronize(backend)
    return field
end

@inline function ka_copy_shifted_complex!(
    out::AbstractMatrix{<:Complex},
    field::AbstractMatrix{<:Complex},
    r0::Int,
    c0::Int,
)
    ny, nx = size(field)
    oy, ox = size(out)
    backend = AK.get_backend(out)
    _ka_copy_shifted_complex_kernel!(backend, (16, 16))(out, field, r0, c0, ny ÷ 2, nx ÷ 2, ny, nx, oy, ox; ndrange=(oy, ox))
    AK.synchronize(backend)
    return out
end

@inline function ka_copy_shifted_intensity!(
    out::AbstractMatrix,
    field::AbstractMatrix{<:Complex},
    r0::Int,
    c0::Int,
)
    ny, nx = size(field)
    oy, ox = size(out)
    backend = AK.get_backend(out)
    _ka_copy_shifted_intensity_kernel!(backend, (16, 16))(out, field, r0, c0, ny ÷ 2, nx ÷ 2, ny, nx, oy, ox; ndrange=(oy, ox))
    AK.synchronize(backend)
    return out
end
