"""Add wavefront phase error in meters, matching upstream PROPER semantics."""
function prop_add_phase(wf::WaveFront, phase_error::Real)
    wf.field .*= cis.((2pi / wf.wavelength_m) * float(phase_error))
    return wf
end

"""Add wavefront phase-error map in meters at current sampling."""
function prop_add_phase(wf::WaveFront, phase_error::AbstractMatrix)
    size(phase_error) == size(wf.field) || throw(ArgumentError("phase size must match wavefront"))
    ph = shift_center_for_wavefront!(wf, phase_error; inverse=true)
    @inbounds wf.field .*= cis.((2pi / wf.wavelength_m) .* ph)
    return wf
end

# Private preparation policy; ordinary prop_add_phase remains the parity oracle.
# Each backend qualifies this arithmetic independently of propagation caching.
_phase_f32_enabled(field) = false
_prepare_phase(wf::WaveFront) = nothing
function _prepare_phase(wf::WaveFront{Float32})
    return _phase_f32_enabled(wf.field) ?
        _PreparedPhase(wf.wavelength_m, 2pi / wf.wavelength_m) : nothing
end

struct _PreparedPhase
    wavelength_m::Float32
    coefficient::Float64
end

@inline function _phase_value_f32(value::ComplexF32, opd::Float32, coefficient::Float64)
    # Retain the original angle, including its Float64 product. Splitting it
    # avoids rounding away the residual before the Float32 trig evaluation.
    angle = coefficient * Float64(opd)
    high = Float32(angle)
    low = Float32(angle - Float64(high))
    magnitude = max(abs(real(value)), abs(imag(value)))
    if isfinite(high) && abs(high) <= 1024f0 && isfinite(low) && abs(low) <= Float32(2^-13) &&
            isfinite(magnitude) && (iszero(magnitude) ||
            floatmin(Float32) <= magnitude <= floatmax(Float32)/4)
        # exp(i*low) ≈ 1+i*low: the leading truncation error is at most 2^-27.
        # Ordinary Float32 trig/product roundoff is validated separately.
        sine, cosine = sincos(high)
        return value * ComplexF32(fma(-low, sine, cosine), fma(low, cosine, sine))
    end
    # Preserve the complete original promoted product outside the fast domain.
    return ComplexF32(value * cis(angle))
end

@kernel function _apply_phase_kernel!(field, phase, coefficient)
    index = @index(Global, Linear)
    if index <= length(field)
        @inbounds field[index] = _phase_value_f32(field[index], phase[index], coefficient)
    end
end

_apply_phase!(wf, phase_error, ::Nothing) = prop_add_phase(wf, phase_error)
function _apply_phase!(wf::WaveFront{Float32}, phase_error::AbstractMatrix{Float32}, prepared::_PreparedPhase)
    size(phase_error) == size(wf.field) || throw(ArgumentError("phase size must match wavefront"))
    wf.wavelength_m == prepared.wavelength_m || throw(ArgumentError("prepared phase wavelength changed"))
    Base.require_one_based_indexing(wf.field, phase_error)
    phase = shift_center_for_wavefront!(wf, phase_error; inverse=true)
    _apply_phase_kernel!(get_backend(wf.field), 256)(wf.field, phase, prepared.coefficient;
        ndrange=length(wf.field))
    return wf
end
