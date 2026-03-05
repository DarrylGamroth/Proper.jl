"""Waist-to-spherical propagation (inside Rayleigh to outside)."""
@inline function _prop_wts_fft_step!(
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

@inline function _prop_wts_fft_step!(
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

function prop_wts(wf::WaveFront, dz::Real, ctx::RunContext, ws::FFTWorkspace)
    wf.reference_surface = SPHERICAL
    iszero(dz) && return wf

    d = float(dz)
    n = size(wf.field, 1)

    wf.z_m += d
    prop_qphase(wf, d)

    _prop_wts_fft_step!(fft_style(ctx), wf, d, n, ctx, ws)

    wf.sampling_m = wf.wavelength_m * abs(d) / (n * wf.sampling_m)
    return wf
end

@inline function prop_wts(wf::WaveFront, dz::Real, ctx::RunContext)
    return prop_wts(wf, dz, ctx, fft_workspace(ctx))
end

@inline function prop_wts(wf::WaveFront, dz::Real)
    return prop_wts(wf, dz, RunContext(wf))
end
