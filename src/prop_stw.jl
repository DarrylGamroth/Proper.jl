"""Spherical-to-waist propagation (outside Rayleigh to inside)."""
@inline function _prop_stw_fft_step!(
    ::FFTWStyle,
    wf::WaveFront{T,<:StridedMatrix{Complex{T}}},
    d::T,
    n::Int,
    ctx::RunContext,
    ws::FFTWorkspace{T},
) where {T<:AbstractFloat}
    ny, nx = size(wf.field)
    f = ensure_fft_scratch!(ws, ny, nx)
    pfft, pbfft = ensure_fft_plans!(ws, ny, nx)
    copyto!(f, wf.field)

    if d >= zero(T)
        LinearAlgebra.mul!(f, pfft, f)
    else
        LinearAlgebra.mul!(f, pbfft, f)
    end

    # Match legacy normalization: both branches scale by 1/n.
    f ./= n
    copyto!(wf.field, f)
    return wf
end

@inline function _prop_stw_fft_step!(
    ::CUFFTStyle,
    wf::WaveFront,
    d::Real,
    n::Int,
    ctx::RunContext,
    ws::FFTWorkspace,
)
    ny, nx = size(wf.field)
    f = ensure_fft_scratch!(ws, ny, nx)
    copyto!(f, wf.field)
    if d >= 0
        fft_forward!(f, ctx)
    else
        fft_backward!(f, ctx)
    end
    f ./= n
    copyto!(wf.field, f)
    return wf
end

@inline function _prop_stw_fft_step!(
    ::FFTStyle,
    wf::WaveFront,
    d::Real,
    n::Int,
    ctx::RunContext,
    ws::FFTWorkspace,
)
    if d >= 0
        f = fft_forward(wf.field, ctx) ./ length(wf.field)
        wf.field .= f .* n
    else
        invf = fft_inverse(wf.field, ctx) .* length(wf.field)
        wf.field .= invf ./ n
    end
    return wf
end

function prop_stw(wf::WaveFront, dz::Real, ctx::RunContext, ws::FFTWorkspace)
    if wf.reference_surface !== SPHERICAL
        return prop_ptp(wf, dz, ctx)
    end

    d = iszero(dz) ? (wf.z_w0_m - wf.z_m) : float(dz)
    n = size(wf.field, 1)

    wf.z_m += d
    wf.sampling_m = wf.wavelength_m * abs(d) / (n * wf.sampling_m)

    _prop_stw_fft_step!(fft_style(ctx), wf, d, n, ctx, ws)

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
