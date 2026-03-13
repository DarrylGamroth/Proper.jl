"""Box-average downsample by integer factor."""
function _prop_pixellate!(
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

function _prop_pixellate!(
    ::SamplingKAExecStyle,
    out::AbstractMatrix,
    img::AbstractMatrix,
    f::Int,
)
    return ka_pixellate!(out, img, f)
end

function _prop_pixellate!(
    out::AbstractMatrix,
    img::AbstractMatrix,
    f::Int,
)
    return _prop_pixellate!(sampling_exec_style(typeof(out), size(out, 1), size(out, 2)), out, img, f)
end

function prop_pixellate!(
    out::AbstractMatrix,
    img::AbstractMatrix,
    factor::Integer,
)
    f = Int(factor)
    f > 0 || throw(ArgumentError("factor must be positive"))
    ny, nx = size(img)
    expected = (ny ÷ f, nx ÷ f)
    size(out) == expected || throw(ArgumentError("output size must be $(expected)"))
    return _prop_pixellate!(out, img, f)
end

function prop_pixellate(img::AbstractMatrix, factor::Integer)
    f = Int(factor)
    f > 0 || throw(ArgumentError("factor must be positive"))
    ny, nx = size(img)
    out = similar(img, ny ÷ f, nx ÷ f)
    return _prop_pixellate!(out, img, f)
end
