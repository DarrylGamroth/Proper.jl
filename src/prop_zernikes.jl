"""Compute simple Zernike-like basis placeholder on wavefront grid."""
function prop_zernikes(wf::WaveFront, nterms::Integer)
    n = Int(nterms)
    ny, nx = size(wf.field)
    out = Array{Float64}(undef, ny, nx, n)
    r = prop_radius(wf)
    rnorm = r ./ maximum(r)
    for k in 1:n
        out[:, :, k] = rnorm .^ (k - 1)
    end
    return out
end
