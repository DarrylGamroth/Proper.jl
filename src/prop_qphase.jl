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
    rsqr = fft_order_rsqr_map(ny, nx, wf.sampling_m)
    wf.field .*= cis.((pi / (wf.wavelength_m * float(c))) .* rsqr)
    return wf
end

"""Apply quadratic phase for curvature c (meters)."""
function prop_qphase(wf::WaveFront, c::Real, ctx::RunContext)
    if wf.field isa StridedMatrix
        return _prop_qphase_strided!(wf, c)
    end
    _ = ctx
    return _prop_qphase_generic!(wf, c)
end

"""Apply quadratic phase for curvature c (meters)."""
function prop_qphase(wf::WaveFront, c::Real)
    return prop_qphase(wf, c, RunContext(typeof(wf.field)))
end
