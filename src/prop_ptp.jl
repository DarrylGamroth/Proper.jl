"""Planar-to-planar Fresnel propagation while keeping planar reference."""
function prop_ptp(wf::WaveFront, dz::Real, ctx::RunContext, ws::FFTWorkspace)
    abs(dz) < 1e-12 && return wf
    wf.reference_surface === PLANAR || throw(ArgumentError("PTP: input reference surface must be PLANAR"))

    n = size(wf.field, 1)
    λ = wf.wavelength_m
    dx = wf.sampling_m

    wf.z_m += float(dz)

    f = fft_forward(wf.field, ctx) ./ length(wf.field)
    f .*= n

    rho2 = ensure_rho2_map!(ws, size(wf.field, 1), size(wf.field, 2), dx)
    f .*= cis.((-pi * λ * float(dz)) .* rho2)

    wf.field .= fft_inverse(f, ctx) .* length(wf.field)
    wf.field ./= n

    return wf
end

@inline function prop_ptp(wf::WaveFront, dz::Real, ctx::RunContext)
    return prop_ptp(wf, dz, ctx, fft_workspace(ctx))
end

@inline function prop_ptp(wf::WaveFront, dz::Real)
    return prop_ptp(wf, dz, RunContext(typeof(wf.field)))
end
