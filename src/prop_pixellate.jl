"""Box-average downsample by integer factor."""
function _prop_pixellate_factor!(
    ::SamplingLoopExecStyle,
    out::AbstractMatrix,
    img::AbstractMatrix,
    f::Int,
)
    ny2, nx2 = size(out)
    scale = inv(float(typeof(real(zero(eltype(out)))))(f * f))
    @inbounds for j in 1:nx2
        for i in 1:ny2
            ys = (i - 1) * f + 1
            xs = (j - 1) * f + 1
            acc = zero(eltype(out))
            for xoff in 0:(f - 1)
                for yoff in 0:(f - 1)
                    acc += img[ys + yoff, xs + xoff]
                end
            end
            out[i, j] = acc * scale
        end
    end
    return out
end

function _prop_pixellate_factor!(
    ::SamplingKAExecStyle,
    out::AbstractMatrix,
    img::AbstractMatrix,
    f::Int,
)
    return ka_pixellate!(out, img, f)
end

function _prop_pixellate_factor!(
    out::AbstractMatrix,
    img::AbstractMatrix,
    f::Int,
)
    return _prop_pixellate_factor!(sampling_exec_style(typeof(out), size(out, 1), size(out, 2)), out, img, f)
end

@inline function _pixellate_axis(n::Integer, center::Integer, mag::T) where {T<:AbstractFloat}
    axis = Vector{T}(undef, n)
    denom = T(center - 1) * T(2) * mag
    shifted = FFTW.ifftshift(collect(1:n) .- center)
    @inbounds for i in eachindex(axis)
        axis[i] = T(shifted[i]) / denom
    end
    return axis
end

function _pixel_mtf(img::AbstractMatrix, sampling_in::Real, sampling_out::Real)
    ny, nx = size(img)
    cy = (ny ÷ 2) + 1
    cx = (nx ÷ 2) + 1
    T = float(promote_type(typeof(sampling_in), typeof(sampling_out), typeof(real(zero(eltype(img))))))
    mag = T(sampling_in / sampling_out)
    vy = _pixellate_axis(ny, cy, mag)
    vx = _pixellate_axis(nx, cx, mag)
    pmtf = Matrix{T}(undef, ny, nx)
    @inbounds for j in 1:nx
        sx = sinc(vx[j])
        for i in 1:ny
            pmtf[i, j] = sinc(vy[i]) * sx
        end
    end
    return pmtf, mag
end

function _prop_pixellate_resample(img::AbstractMatrix, sampling_in::Real, sampling_out::Real, n_out::Integer)
    pmtf, mag = _pixel_mtf(img, sampling_in, sampling_out)
    shifted = FFTW.ifftshift(img)
    amtf = pmtf .* fft(shifted)
    cimc = FFTW.fftshift(abs.(ifft(amtf)) ./ (mag * mag))
    out_n = n_out > 0 ? Int(n_out) : floor(Int, size(img, 2) * mag)
    return prop_magnify(cimc, mag, out_n)
end

function _prop_pixellate_factor!(
    out::AbstractMatrix,
    img::AbstractMatrix,
    factor::Integer,
)
    f = Int(factor)
    f > 0 || throw(ArgumentError("factor must be positive"))
    ny, nx = size(img)
    expected = (ny ÷ f, nx ÷ f)
    size(out) == expected || throw(ArgumentError("output size must be $(expected)"))
    return _prop_pixellate_factor!(out, img, f)
end

function _prop_pixellate_factor(
    img::AbstractMatrix,
    factor::Integer,
)
    f = Int(factor)
    f > 0 || throw(ArgumentError("factor must be positive"))
    ny, nx = size(img)
    out = similar(img, ny ÷ f, nx ÷ f)
    return _prop_pixellate_factor!(out, img, f)
end

"""
    prop_pixellate(img, sampling_in, sampling_out, n_out=0)

Integrate a sampled PSF onto square detector pixels.

This matches the upstream PROPER public API. The input PSF is convolved with
the transfer function of an ideal square pixel, then resampled to detector
pixel spacing. `sampling_in` and `sampling_out` are in meters per pixel. If
`n_out == 0`, the output size follows the magnification implied by the sampling
ratio.
"""
function prop_pixellate(
    img::AbstractMatrix,
    sampling_in::Real,
    sampling_out::Real,
    n_out::Integer=0,
)
    sampling_in > 0 || throw(ArgumentError("sampling_in must be positive"))
    sampling_out > 0 || throw(ArgumentError("sampling_out must be positive"))
    return _prop_pixellate_resample(img, sampling_in, sampling_out, n_out)
end

"""
    prop_pixellate!(out, img, sampling_in, sampling_out)

Integrate a sampled PSF onto detector pixels and write the result into `out`.

`out` must match the size implied by the requested detector sampling and chosen
output dimensions.
"""
function prop_pixellate!(
    out::AbstractMatrix,
    img::AbstractMatrix,
    sampling_in::Real,
    sampling_out::Real,
)
    sampling_in > 0 || throw(ArgumentError("sampling_in must be positive"))
    sampling_out > 0 || throw(ArgumentError("sampling_out must be positive"))
    tmp = _prop_pixellate_resample(img, sampling_in, sampling_out, size(out, 2))
    size(tmp) == size(out) || throw(ArgumentError("output size must be $(size(tmp))"))
    copyto!(out, tmp)
    return out
end
