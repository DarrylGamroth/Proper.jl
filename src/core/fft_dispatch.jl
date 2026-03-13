@inline function fft_forward(a::AbstractArray, ctx::RunContext)
    return _fft_forward(fft_style(ctx), a)
end

@inline function fft_forward!(a::AbstractArray, ctx::RunContext)
    return _fft_forward!(fft_style(ctx), a)
end

@inline function fft_inverse(a::AbstractArray, ctx::RunContext)
    return _fft_inverse(fft_style(ctx), a)
end

@inline function fft_backward!(a::AbstractArray, ctx::RunContext)
    return _fft_backward!(fft_style(ctx), a)
end

@inline _fft_forward(::FFTWStyle, a::AbstractArray) = fft(a)
@inline _fft_forward(::CUFFTStyle, a::AbstractArray) = fft(a)
@inline _fft_forward(::FFTStyle, a::AbstractArray) = fft(a)

@inline _fft_forward!(::FFTWStyle, a::AbstractArray) = fft!(a)
@inline _fft_forward!(::CUFFTStyle, a::AbstractArray) = fft!(a)
@inline _fft_forward!(::FFTStyle, a::AbstractArray) = fft!(a)

@inline _fft_inverse(::FFTWStyle, a::AbstractArray) = ifft(a)
@inline _fft_inverse(::CUFFTStyle, a::AbstractArray) = ifft(a)
@inline _fft_inverse(::FFTStyle, a::AbstractArray) = ifft(a)

@inline _fft_backward!(::FFTWStyle, a::AbstractArray) = bfft!(a)
@inline _fft_backward!(::CUFFTStyle, a::AbstractArray) = bfft!(a)
@inline _fft_backward!(::FFTStyle, a::AbstractArray) = bfft!(a)
