abstract type PTPExecStyle end
struct PTPPlannedExecStyle <: PTPExecStyle end
struct PTPKAExecStyle <: PTPExecStyle end
struct PTPGenericExecStyle <: PTPExecStyle end

@inline ptp_exec_style(::StridedLayout, ::CPUBackend, ::FFTWStyle) = PTPPlannedExecStyle()
@inline ptp_exec_style(::ArrayLayoutStyle, ::CUDABackend, ::CUFFTStyle) = PTPKAExecStyle()
@inline ptp_exec_style(::ArrayLayoutStyle, ::BackendStyle, ::FFTStyle) = PTPGenericExecStyle()

function _prop_ptp_fft!(
    ::PTPPlannedExecStyle,
    wf::WaveFront,
    ctx::RunContext,
    ws::FFTWorkspace,
    n::Int,
    ny::Int,
    nx::Int,
    dx,
    kphase,
)
    f = prepare_fft_field!(ws, wf.field, ny, nx)
    pfft, pbfft = ensure_fft_plans!(ws, ny, nx, fft_planning_style(ctx))
    LinearAlgebra.mul!(f, pfft, f)
    T = typeof(float(dx))
    _apply_fft_quadratic_phase_strided!(f, T(kphase), T(dx), qphase_workspace(ctx), inv(T(n)))

    LinearAlgebra.mul!(f, pbfft, f)
    f ./= n
    wf.field = f
    return wf
end

function _prop_ptp_fft!(
    ::PTPKAExecStyle,
    wf::WaveFront,
    ctx::RunContext,
    ws::FFTWorkspace,
    n::Int,
    ny::Int,
    nx::Int,
    dx,
    kphase,
)
    f = prepare_fft_field!(ws, wf.field, ny, nx)
    pfft, pbfft = ensure_fft_plans!(ws, ny, nx, fft_planning_style(ctx))
    LinearAlgebra.mul!(f, pfft, f)
    f ./= n
    ka_apply_frequency_phase!(f, kphase, dx)
    LinearAlgebra.mul!(f, pbfft, f)
    f ./= n
    wf.field = f
    return wf
end

function _prop_ptp_fft!(
    ::PTPGenericExecStyle,
    wf::WaveFront,
    ctx::RunContext,
    ws::FFTWorkspace,
    n::Int,
    ny::Int,
    nx::Int,
    dx,
    kphase,
)
    f = fft_forward(wf.field, ctx) ./ length(wf.field)
    f .*= n

    rho2 = backend_adapt(f, ensure_rho2_map!(ws, ny, nx, dx))
    f .*= cis.(kphase .* rho2)

    wf.field .= fft_inverse(f, ctx) .* length(wf.field)
    wf.field ./= n
    return wf
end

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
    sty = ptp_exec_style(array_layout_style(typeof(wf.field)), ctx.backend, fft_style(ctx))
    return _prop_ptp_fft!(sty, wf, ctx, ws, n, ny, nx, dx, kphase)
end

@inline function prop_ptp(wf::WaveFront, dz::Real, ctx::RunContext)
    return prop_ptp(wf, dz, ctx, fft_workspace(ctx))
end

@inline function prop_ptp(wf::WaveFront, dz::Real)
    return prop_ptp(wf, dz, RunContext(wf))
end
