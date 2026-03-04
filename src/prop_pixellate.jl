"""Box-average downsample by integer factor."""
function prop_pixellate(img::AbstractMatrix, factor::Integer)
    f = Int(factor)
    f > 0 || throw(ArgumentError("factor must be positive"))
    ny, nx = size(img)
    ny2 = ny ÷ f
    nx2 = nx ÷ f
    out = similar(img, ny2, nx2)
    @inbounds for j in 1:nx2
        for i in 1:ny2
            ys = (i - 1) * f + 1
            xs = (j - 1) * f + 1
            out[i, j] = sum(@view img[ys:(ys + f - 1), xs:(xs + f - 1)]) / (f * f)
        end
    end
    return out
end
