"""Planar-to-planar Fresnel propagation while keeping planar reference."""
function prop_ptp(wf::WaveFront, dz::Real, ctx::RunContext, ws::FFTWorkspace)
    abs(dz) < 1e-12 && return wf
    wf.reference_surface === PLANAR || throw(ArgumentError("PTP: input reference surface must be PLANAR"))

    n = size(wf.field, 1)
    λ = wf.wavelength_m
    dx = wf.sampling_m

    wf.z_m += float(dz)

    ny, nx = size(wf.field)
    kphase = -pi * λ * float(dz)

    if wf.field isa StridedMatrix && fft_style(ctx) isa FFTWStyle
        f = ensure_fft_scratch!(ws, ny, nx)
        pfft, pbfft = ensure_fft_plans!(ws, ny, nx)
        copyto!(f, wf.field)
        LinearAlgebra.mul!(f, pfft, f)
        f ./= n

        rho2 = ensure_rho2_map!(ws, ny, nx, dx)
        @inbounds @simd for idx in eachindex(f, rho2)
            f[idx] *= cis(kphase * rho2[idx])
        end

        LinearAlgebra.mul!(f, pbfft, f)
        f ./= n
        copyto!(wf.field, f)
    else
        f = fft_forward(wf.field, ctx) ./ length(wf.field)
        f .*= n

        rho2 = backend_adapt(f, ensure_rho2_map!(ws, ny, nx, dx))
        f .*= cis.(kphase .* rho2)

        wf.field .= fft_inverse(f, ctx) .* length(wf.field)
        wf.field ./= n
    end

    return wf
end

@inline function prop_ptp(wf::WaveFront, dz::Real, ctx::RunContext)
    return prop_ptp(wf, dz, ctx, fft_workspace(ctx))
end

@inline function prop_ptp(wf::WaveFront, dz::Real)
    return prop_ptp(wf, dz, RunContext(wf))
end
