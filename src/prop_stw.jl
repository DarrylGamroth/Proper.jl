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
    f = prepare_fft_field!(ws, wf.field, ny, nx)
    pfft, pbfft = ensure_fft_plans!(ws, ny, nx, fft_planning_style(ctx))
    padded = _padded_propagation_ffts!(ws, pfft, pbfft)

    if d >= zero(T)
        _propagation_fft!(f, pfft, padded)
    else
        _propagation_fft!(f, pbfft, padded)
    end

    # Match legacy normalization: both branches scale by 1/n.
    rmul!(f, inv(T(n)))
    wf.field = f
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

@inline function _prop_stw_fft_step!(
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

@inline function _prop_stw_fft_step!(
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

@inline function _prop_stw!(wf::OpticalSurfaceState, dz::Real, ctx, ws)
    n = _require_square_propagation_grid(wf, :prop_stw)
    if wf.reference_surface !== SPHERICAL
        return _execute_propagation!(prop_ptp, wf, dz, ctx)
    end

    d = iszero(dz) ? (wf.z_w0_m - wf.z_m) : float(dz)

    wf.z_m += d
    wf.sampling_m = wf.wavelength_m * abs(d) / (n * wf.sampling_m)

    _stw_field!(wf, d, ctx, ws, n)
    wf.reference_surface = PLANAR
    return wf
end

@inline _stw_field!(wf::WaveFrontGeometry, d, ::Nothing, ::Nothing, n) = _check_spatial_phase(wf, d)

@inline function _stw_field!(wf::WaveFront, d, ctx::RunContext, ws::FFTWorkspace, n)
    _prop_stw_fft_step!(fft_style(ctx), wf, d, n, ctx, ws)
    prop_qphase(wf, d, ctx)
    return apply_carrier_phase!(wf, d, ctx)
end

@inline prop_stw(wf::WaveFront, dz::Real, ctx::RunContext, ws::FFTWorkspace) =
    _prop_stw!(wf, dz, ctx, ws)

@inline prop_stw(wf::WaveFrontGeometry, dz::Real=0.0) = _prop_stw!(wf, dz, nothing, nothing)

"""Propagate from a spherical-reference regime back to a waist/planar regime."""
@inline function prop_stw(wf::WaveFront, dz::Real, ctx::RunContext)
    return prop_stw(wf, dz, ctx, fft_workspace(ctx))
end

"""Propagate from a spherical-reference regime back to a waist/planar regime."""
@inline function prop_stw(wf::WaveFront, dz::Real=0.0)
    return prop_stw(wf, dz, resolve_run_context(wf))
end
