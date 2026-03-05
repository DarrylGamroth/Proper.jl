"""Spherical-to-waist propagation (outside Rayleigh to inside)."""
function prop_stw(wf::WaveFront, dz::Real, ctx::RunContext)
    if wf.reference_surface !== SPHERICAL
        return prop_ptp(wf, dz, ctx)
    end

    d = iszero(dz) ? (wf.z_w0_m - wf.z_m) : float(dz)
    n = size(wf.field, 1)

    wf.z_m += d
    wf.sampling_m = wf.wavelength_m * abs(d) / (n * wf.sampling_m)

    if d >= 0
        f = fft_forward(wf.field, ctx) ./ length(wf.field)
        wf.field .= f .* n
    else
        invf = fft_inverse(wf.field, ctx) .* length(wf.field)
        wf.field .= invf ./ n
    end

    prop_qphase(wf, d)
    wf.reference_surface = PLANAR
    return wf
end

@inline function prop_stw(wf::WaveFront, dz::Real=0.0)
    return prop_stw(wf, dz, RunContext(typeof(wf.field)))
end
