"""Waist-to-spherical propagation (inside Rayleigh to outside)."""
function prop_wts(wf::WaveFront, dz::Real, ctx::RunContext)
    wf.reference_surface = SPHERICAL
    iszero(dz) && return wf

    d = float(dz)
    n = size(wf.field, 1)

    wf.z_m += d
    prop_qphase(wf, d)

    if d >= 0
        f = fft_forward(wf.field, ctx) ./ length(wf.field)
        wf.field .= f .* n
    else
        invf = fft_inverse(wf.field, ctx) .* length(wf.field)
        wf.field .= invf ./ n
    end

    wf.sampling_m = wf.wavelength_m * abs(d) / (n * wf.sampling_m)
    return wf
end

@inline function prop_wts(wf::WaveFront, dz::Real)
    return prop_wts(wf, dz, RunContext(typeof(wf.field)))
end
