function ffts!(wavefront::AbstractMatrix{<:Complex}, direction::Integer)
    n = size(wavefront, 1)
    shifted = circshift(wavefront, (-n ÷ 2, -n ÷ 2))
    transformed = if direction == -1
        fft(shifted) ./ length(shifted)
    else
        ifft(shifted) .* length(shifted)
    end
    copyto!(wavefront, circshift(transformed, (n ÷ 2, n ÷ 2)))
    return wavefront
end

function ffts(wavefront::AbstractMatrix, direction::Integer)
    arr = wavefront isa AbstractMatrix{<:Complex} ? copy(wavefront) : ComplexF64.(wavefront)
    return ffts!(arr, direction)
end
