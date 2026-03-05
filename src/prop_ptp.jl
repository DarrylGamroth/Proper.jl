"""Planar-to-planar Fresnel propagation while keeping planar reference."""
function prop_ptp(wf::WaveFront, dz::Real)
    abs(dz) < 1e-12 && return wf
    wf.reference_surface == :PLANAR || throw(ArgumentError("PTP: input reference surface must be PLANAR"))

    n = size(wf.field, 1)
    λ = wf.wavelength_m
    dx = wf.sampling_m

    wf.z_m += float(dz)

    f = fft(wf.field) ./ length(wf.field)
    f .*= n

    rho2 = fft_order_rho2_map(size(wf.field, 1), size(wf.field, 2), dx)
    f .*= cis.((-pi * λ * float(dz)) .* rho2)

    wf.field .= ifft(f) .* length(wf.field)
    wf.field ./= n

    return wf
end
