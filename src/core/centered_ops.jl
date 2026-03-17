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
