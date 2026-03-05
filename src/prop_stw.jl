"""Spherical-to-waist propagation (outside Rayleigh to inside)."""
function prop_stw(wf::WaveFront, dz::Real, ctx::RunContext, ws::FFTWorkspace)
    if wf.reference_surface !== SPHERICAL
        return prop_ptp(wf, dz, ctx)
    end

    d = iszero(dz) ? (wf.z_w0_m - wf.z_m) : float(dz)
    ny, nx = size(wf.field)
    n = ny

    wf.z_m += d
    wf.sampling_m = wf.wavelength_m * abs(d) / (n * wf.sampling_m)

    if wf.field isa StridedMatrix && fft_style(ctx) isa FFTWStyle
        f = ensure_fft_scratch!(ws, ny, nx)
        pfft, pbfft = ensure_fft_plans!(ws, ny, nx)
        copyto!(f, wf.field)

        if d >= 0
            LinearAlgebra.mul!(f, pfft, f)
        else
            LinearAlgebra.mul!(f, pbfft, f)
        end

        # Match legacy normalization: both branches scale by 1/n.
        f ./= n
        copyto!(wf.field, f)
    else
        if d >= 0
            f = fft_forward(wf.field, ctx) ./ length(wf.field)
            wf.field .= f .* n
        else
            invf = fft_inverse(wf.field, ctx) .* length(wf.field)
            wf.field .= invf ./ n
        end
    end

    prop_qphase(wf, d)
    wf.reference_surface = PLANAR
    return wf
end

@inline function prop_stw(wf::WaveFront, dz::Real, ctx::RunContext)
    return prop_stw(wf, dz, ctx, fft_workspace(ctx))
end

@inline function prop_stw(wf::WaveFront, dz::Real=0.0)
    return prop_stw(wf, dz, RunContext(wf))
end
