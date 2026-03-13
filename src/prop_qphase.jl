@inline function _shifted_index_0based(p::Int, n::Int)
    c = n ÷ 2
    cut = n - c - 1
    return p <= cut ? p : (p - n)
end

@inline function _prop_qphase_strided!(wf::WaveFront, c::Real)
    c == 0 && return wf
    ny, nx = size(wf.field)
    dx = wf.sampling_m
    k = pi / (wf.wavelength_m * float(c))

    @inbounds for j in 1:nx
        x0 = _shifted_index_0based(j - 1, nx)
        x2 = (x0 * dx)^2
        for i in 1:ny
            y0 = _shifted_index_0based(i - 1, ny)
            rsqr = x2 + (y0 * dx)^2
            wf.field[i, j] *= cis(k * rsqr)
        end
    end

    return wf
end

@inline function _prop_qphase_generic!(wf::WaveFront, c::Real)
    c == 0 && return wf
    ny, nx = size(wf.field)
    rsqr = backend_adapt(wf.field, fft_order_rsqr_map(ny, nx, wf.sampling_m))
    wf.field .*= cis.((pi / (wf.wavelength_m * float(c))) .* rsqr)
    return wf
end

abstract type QPhaseExecStyle end
struct QPhaseLoopExecStyle <: QPhaseExecStyle end
struct QPhaseBroadcastExecStyle <: QPhaseExecStyle end

@inline qphase_exec_style(::StridedLayout, ::CPUBackend) = QPhaseLoopExecStyle()
@inline qphase_exec_style(::ArrayLayoutStyle, ::BackendStyle) = QPhaseBroadcastExecStyle()

@inline _prop_qphase!(::QPhaseLoopExecStyle, wf::WaveFront, c::Real) = _prop_qphase_strided!(wf, c)
@inline _prop_qphase!(::QPhaseBroadcastExecStyle, wf::WaveFront, c::Real) = _prop_qphase_generic!(wf, c)

"""Apply quadratic phase for curvature c (meters)."""
function prop_qphase(wf::WaveFront, c::Real, ctx::RunContext)
    sty = qphase_exec_style(array_layout_style(typeof(wf.field)), ctx.backend)
    return _prop_qphase!(sty, wf, c)
end

"""Apply quadratic phase for curvature c (meters)."""
function prop_qphase(wf::WaveFront, c::Real)
    return prop_qphase(wf, c, RunContext(wf))
end
