@inline function _shifted_index_0based(p::Int, n::Int)
    c = n ÷ 2
    cut = n - c - 1
    return p <= cut ? p : (p - n)
end

@inline function _apply_separable_phase!(field::AbstractMatrix{Complex{T}}, xphase::AbstractVector{Complex{T}}, yphase::AbstractVector{Complex{T}}, scale::T=one(T)) where {T<:AbstractFloat}
    ny, nx = size(field)
    @inbounds for j in 1:nx
        xj = scale * xphase[j]
        for i in 1:ny
            field[i, j] *= xj * yphase[i]
        end
    end
    return field
end

@inline function _fill_fft_order_axis_phase!(phase::AbstractVector{Complex{T}}, n::Int, spacing::T, k::T) where {T<:AbstractFloat}
    @inbounds for p in 1:n
        x0 = _shifted_index_0based(p - 1, n) * spacing
        phase[p] = cis(k * x0 * x0)
    end
    return phase
end

@inline function _prop_qphase_strided!(wf::WaveFront, c::Real, ws::QPhaseWorkspace, scale::Real=1)
    T = float(real(eltype(wf.field)))
    c == 0 && return isone(scale) ? wf : (rmul!(wf.field, T(scale)); wf)
    ny, nx = size(wf.field)
    dx = wf.sampling_m
    k = pi / (wf.wavelength_m * float(c))
    xphase, yphase = ensure_qphase_vectors!(ws, nx, ny)

    _fill_fft_order_axis_phase!(xphase, nx, dx, k)
    _fill_fft_order_axis_phase!(yphase, ny, dx, k)
    _apply_separable_phase!(wf.field, xphase, yphase, T(scale))

    return wf
end

@inline function _apply_fft_quadratic_phase_strided!(
    field::AbstractMatrix{Complex{T}},
    kphase::T,
    dx::T,
    ws::QPhaseWorkspace{T},
    scale::T=one(T),
) where {T<:AbstractFloat}
    ny, nx = size(field)
    inv_dy = inv(T(ny) * dx)
    inv_dx = inv(T(nx) * dx)
    xphase, yphase = ensure_qphase_vectors!(ws, nx, ny)
    _fill_fft_order_axis_phase!(xphase, nx, inv_dx, kphase)
    _fill_fft_order_axis_phase!(yphase, ny, inv_dy, kphase)
    _apply_separable_phase!(field, xphase, yphase, scale)
    return field
end

@inline function _prop_qphase_generic!(wf::WaveFront, c::Real)
    c == 0 && return wf
    ny, nx = size(wf.field)
    rsqr = backend_adapt(wf.field, fft_order_rsqr_map(ny, nx, wf.sampling_m))
    wf.field .*= cis.((pi / (wf.wavelength_m * float(c))) .* rsqr)
    return wf
end

@inline function _prop_qphase_ka!(wf::WaveFront, c::Real)
    c == 0 && return wf
    k = pi / (wf.wavelength_m * float(c))
    ka_apply_qphase!(wf.field, k, wf.sampling_m)
    return wf
end

abstract type QPhaseExecStyle end
struct QPhaseLoopExecStyle <: QPhaseExecStyle end
struct QPhaseKAExecStyle <: QPhaseExecStyle end
struct QPhaseBroadcastExecStyle <: QPhaseExecStyle end

@inline qphase_exec_style(::StridedLayout, ::CPUBackend) = QPhaseLoopExecStyle()
@inline qphase_exec_style(::ArrayLayoutStyle, ::CUDABackend) = QPhaseKAExecStyle()
@inline qphase_exec_style(::ArrayLayoutStyle, ::AMDGPUBackend) = QPhaseKAExecStyle()
@inline qphase_exec_style(::ArrayLayoutStyle, ::BackendStyle) = QPhaseBroadcastExecStyle()

@inline _prop_qphase!(::QPhaseLoopExecStyle, wf::WaveFront, c::Real, ctx::RunContext) = _prop_qphase_strided!(wf, c, qphase_workspace(ctx))
@inline _prop_qphase!(::QPhaseKAExecStyle, wf::WaveFront, c::Real, ctx::RunContext) = _prop_qphase_ka!(wf, c)
@inline _prop_qphase!(::QPhaseBroadcastExecStyle, wf::WaveFront, c::Real, ctx::RunContext) = _prop_qphase_generic!(wf, c)

"""Apply quadratic phase for curvature c (meters)."""
function prop_qphase(wf::WaveFront, c::Real, ctx::RunContext)
    sty = qphase_exec_style(array_layout_style(typeof(wf.field)), ctx.backend)
    return _prop_qphase!(sty, wf, c, ctx)
end

"""Apply quadratic phase for curvature c (meters)."""
function prop_qphase(wf::WaveFront, c::Real)
    return prop_qphase(wf, c, RunContext(wf))
end
