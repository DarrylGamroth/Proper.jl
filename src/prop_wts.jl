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
    f = prepare_fft_field!(ws, wf.field, ny, nx)
    pfft, pbfft = ensure_fft_plans!(ws, ny, nx, fft_planning_style(ctx))

    if d >= zero(T)
        LinearAlgebra.mul!(f, pfft, f)
    else
        LinearAlgebra.mul!(f, pbfft, f)
    end

    # Match legacy normalization: both branches scale by 1/n.
    rmul!(f, inv(T(n)))
    wf.field = f
    return wf
end

@inline function _prop_wts_fft_step!(
    ::CUFFTStyle,
    wf::WaveFront,
    d::Real,
    n::Int,
    ctx::RunContext,
    ws::FFTWorkspace,
)
    ny, nx = size(wf.field)
    f = prepare_fft_field!(ws, wf.field, ny, nx)
    pfft, pbfft = ensure_fft_plans!(ws, ny, nx, fft_planning_style(ctx))
    T = float(real(eltype(f)))
    if d >= 0
        LinearAlgebra.mul!(f, pfft, f)
    else
        LinearAlgebra.mul!(f, pbfft, f)
    end
    ka_scale_field!(f, inv(T(n)))
    wf.field = f
    return wf
end

@inline function _prop_wts_fft_step!(
    ::ROCFFTStyle,
    wf::WaveFront,
    d::Real,
    n::Int,
    ctx::RunContext,
    ws::FFTWorkspace,
)
    ny, nx = size(wf.field)
    f = prepare_fft_field!(ws, wf.field, ny, nx)
    pfft, pbfft = ensure_fft_plans!(ws, ny, nx, fft_planning_style(ctx))
    T = float(real(eltype(f)))
    if d >= 0
        pfft * f
    else
        pbfft * f
    end
    ka_scale_field!(f, inv(T(n)))
    wf.field = f
    return wf
end

@inline function _prop_wts_fft_step!(
    ::FFTStyle,
    wf::WaveFront,
    d::Real,
    n::Int,
    ::RunContext,
    ws::FFTWorkspace,
)
    _ = ws
    if d >= 0
        f = fft(wf.field) ./ length(wf.field)
        wf.field .= f .* n
    else
        invf = ifft(wf.field) .* length(wf.field)
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
    prop_qphase(wf, d, ctx)

    _prop_wts_fft_step!(fft_style(ctx), wf, d, n, ctx, ws)

    wf.sampling_m = wf.wavelength_m * abs(d) / (n * wf.sampling_m)
    return wf
end

"""Propagate from a waist/planar regime to a spherical-reference regime."""
@inline function prop_wts(wf::WaveFront, dz::Real, ctx::RunContext)
    return prop_wts(wf, dz, ctx, fft_workspace(ctx))
end

"""Propagate from a waist/planar regime to a spherical-reference regime."""
@inline function prop_wts(wf::WaveFront, dz::Real)
    return prop_wts(wf, dz, RunContext(wf))
end
