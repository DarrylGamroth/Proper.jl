mutable struct CenteredFFTCache{T<:AbstractFloat,FP,BP}
    scratch::Matrix{Complex{T}}
    shifted::Matrix{Complex{T}}
    forward_plan::FP
    backward_plan::BP
end

function CenteredFFTCache(::Type{T}, n::Integer; flags=FFTW.MEASURE) where {T<:AbstractFloat}
    scratch = Matrix{Complex{T}}(undef, n, n)
    shifted = Matrix{Complex{T}}(undef, n, n)
    return CenteredFFTCache{T,typeof(FFTW.plan_fft!(scratch; flags=flags)),typeof(FFTW.plan_bfft!(scratch; flags=flags))}(
        scratch,
        shifted,
        FFTW.plan_fft!(scratch; flags=flags),
        FFTW.plan_bfft!(scratch; flags=flags),
    )
end

function copy_centered!(out::AbstractMatrix{T}, input_image::AbstractMatrix) where {T}
    input_dim = size(input_image, 2)
    output_dim = size(out, 2)
    if input_dim == output_dim
        out .= T.(input_image)
        return out
    elseif output_dim < input_dim
        x1 = input_dim ÷ 2 - output_dim ÷ 2 + 1
        x2 = x1 + output_dim - 1
        @views out .= T.(input_image[x1:x2, x1:x2])
        return out
    end

    fill!(out, zero(T))
    x1 = output_dim ÷ 2 - input_dim ÷ 2 + 1
    x2 = x1 + input_dim - 1
    @views out[x1:x2, x1:x2] .= T.(input_image)
    return out
end

function half_shift_copy!(out::AbstractMatrix{T}, input::AbstractMatrix) where {T}
    return shift_copy!(out, input, size(input, 1) ÷ 2, size(input, 2) ÷ 2)
end

function copy_centered_to_fft_order!(out::AbstractMatrix{T}, input_image::AbstractMatrix) where {T}
    ny_out, nx_out = size(out)
    ny_in, nx_in = size(input_image)
    oy1 = ny_out ÷ 2 - ny_in ÷ 2 + 1
    ox1 = nx_out ÷ 2 - nx_in ÷ 2 + 1
    oy2 = oy1 + ny_in - 1
    ox2 = ox1 + nx_in - 1

    sy1 = 1
    sx1 = 1
    sy2 = ny_in
    sx2 = nx_in

    if oy1 < 1
        sy1 += 1 - oy1
        oy1 = 1
    end
    if ox1 < 1
        sx1 += 1 - ox1
        ox1 = 1
    end
    if oy2 > ny_out
        sy2 -= oy2 - ny_out
        oy2 = ny_out
    end
    if ox2 > nx_out
        sx2 -= ox2 - nx_out
        ox2 = nx_out
    end

    row_split = ny_out - (ny_out ÷ 2)
    col_split = nx_out - (nx_out ÷ 2)
    @inbounds for j in 1:nx_out
        pj = j <= (nx_out ÷ 2) ? (col_split + j) : (j - (nx_out ÷ 2))
        for i in 1:ny_out
            pi = i <= (ny_out ÷ 2) ? (row_split + i) : (i - (ny_out ÷ 2))
            if oy1 <= pi <= oy2 && ox1 <= pj <= ox2
                out[i, j] = T(input_image[sy1 + (pi - oy1), sx1 + (pj - ox1)])
            else
                out[i, j] = zero(T)
            end
        end
    end
    return out
end

function copy_fft_order_resized!(out::AbstractMatrix{T}, input_image::AbstractMatrix) where {T}
    ny_out, nx_out = size(out)
    ny_in, nx_in = size(input_image)
    oy1 = ny_out ÷ 2 - ny_in ÷ 2 + 1
    ox1 = nx_out ÷ 2 - nx_in ÷ 2 + 1
    oy2 = oy1 + ny_in - 1
    ox2 = ox1 + nx_in - 1

    sy1 = 1
    sx1 = 1
    sy2 = ny_in
    sx2 = nx_in

    if oy1 < 1
        sy1 += 1 - oy1
        oy1 = 1
    end
    if ox1 < 1
        sx1 += 1 - ox1
        ox1 = 1
    end
    if oy2 > ny_out
        sy2 -= oy2 - ny_out
        oy2 = ny_out
    end
    if ox2 > nx_out
        sx2 -= ox2 - nx_out
        ox2 = nx_out
    end

    row_split_out = ny_out - (ny_out ÷ 2)
    col_split_out = nx_out - (nx_out ÷ 2)
    sy_shift = ny_in ÷ 2
    sx_shift = nx_in ÷ 2
    @inbounds for j in 1:nx_out
        pj = j <= (nx_out ÷ 2) ? (col_split_out + j) : (j - (nx_out ÷ 2))
        for i in 1:ny_out
            pi = i <= (ny_out ÷ 2) ? (row_split_out + i) : (i - (ny_out ÷ 2))
            if oy1 <= pi <= oy2 && ox1 <= pj <= ox2
                src_i = sy1 + (pi - oy1)
                src_j = sx1 + (pj - ox1)
                out[i, j] = T(input_image[mod1(src_i + sy_shift, ny_in), mod1(src_j + sx_shift, nx_in)])
            else
                out[i, j] = zero(T)
            end
        end
    end
    return out
end

function shift_copy!(out::AbstractMatrix{T}, input::AbstractMatrix, sy::Integer, sx::Integer) where {T}
    ny, nx = size(input)
    size(out) == (ny, nx) || throw(ArgumentError("shift output size must match input size"))
    0 <= sy <= ny || throw(ArgumentError("row shift must be within array bounds"))
    0 <= sx <= nx || throw(ArgumentError("column shift must be within array bounds"))
    top_h = sy
    left_w = sx
    row_split = ny - sy
    col_split = nx - sx

    @views begin
        if sy > 0 && sx > 0
            out[1:top_h, 1:left_w] .= T.(input[(row_split + 1):ny, (col_split + 1):nx])
        end
        if sy > 0 && sx < nx
            out[1:top_h, (left_w + 1):nx] .= T.(input[(row_split + 1):ny, 1:col_split])
        end
        if sy < ny && sx > 0
            out[(top_h + 1):ny, 1:left_w] .= T.(input[1:row_split, (col_split + 1):nx])
        end
        if sy < ny && sx < nx
            out[(top_h + 1):ny, (left_w + 1):nx] .= T.(input[1:row_split, 1:col_split])
        end
    end
    return out
end

function reverse_shift1!(out::AbstractMatrix{Complex{T}}, input::AbstractMatrix{<:Complex}) where {T<:AbstractFloat}
    ny, nx = size(input)
    size(out) == (ny, nx) || throw(ArgumentError("reverse-shift output size must match input size"))
    @inbounds for j in 1:nx
        src_j = mod1(2 - j, nx)
        for i in 1:ny
            src_i = mod1(2 - i, ny)
            out[i, j] = input[src_i, src_j]
        end
    end
    return out
end

function centered_fft!(field::Matrix{Complex{T}}, cache::CenteredFFTCache{T}, direction::Integer) where {T<:AbstractFloat}
    half_shift_copy!(cache.scratch, field)
    if direction == -1
        cache.forward_plan * cache.scratch
        cache.scratch .*= inv(length(cache.scratch))
    else
        cache.backward_plan * cache.scratch
    end
    half_shift_copy!(field, cache.scratch)
    return field
end
