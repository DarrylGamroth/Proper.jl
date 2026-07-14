using Proper

abstract type AbstractOcculter end
struct GaussianOcculter <: AbstractOcculter end
struct SolidOcculter <: AbstractOcculter end
struct EighthOrderOcculter <: AbstractOcculter end

normalize_occulter(occulter::AbstractOcculter) = occulter

function normalize_occulter(occulter::Symbol)
    occulter in (:gaussian, :GaussianOcculter) && return GaussianOcculter()
    occulter in (:solid, :SolidOcculter) && return SolidOcculter()
    occulter in (:eighth_order, :EighthOrderOcculter) && return EighthOrderOcculter()
    throw(ArgumentError("Unknown occulter selector: $(occulter)"))
end

function normalize_occulter(occulter::AbstractString)
    normalized = lowercase(replace(occulter, '-' => '_', ' ' => '_'))
    normalized == "gaussian" && return GaussianOcculter()
    normalized == "solid" && return SolidOcculter()
    normalized in ("8th_order", "eighth_order") && return EighthOrderOcculter()
    throw(ArgumentError("Unknown occulter selector: $(occulter)"))
end

"""Caller-owned, centered amplitude snapshots from the coronagraph example."""
struct CoronagraphDiagnostics{T<:Real,A<:AbstractMatrix{T}}
    after_occulter_amplitude::A
    before_lyot_amplitude::A

    function CoronagraphDiagnostics(
        after_occulter_amplitude::A,
        before_lyot_amplitude::A,
    ) where {T<:Real,A<:AbstractMatrix{T}}
        size(after_occulter_amplitude) == size(before_lyot_amplitude) ||
            throw(ArgumentError("coronagraph diagnostic buffers must have matching sizes"))
        return new{T,A}(after_occulter_amplitude, before_lyot_amplitude)
    end
end

function CoronagraphDiagnostics(wfo::WaveFront{T}) where {T}
    after_occulter_amplitude = similar(wfo.field, T)
    before_lyot_amplitude = similar(wfo.field, T)
    return CoronagraphDiagnostics(after_occulter_amplitude, before_lyot_amplitude)
end

@inline _capture_after_occulter!(::Nothing, ::WaveFront) = nothing
@inline _capture_before_lyot!(::Nothing, ::WaveFront) = nothing

function _check_coronagraph_diagnostic_buffer(buffer::AbstractMatrix, wfo::WaveFront{T}, label::AbstractString) where {T}
    size(buffer) == size(wfo.field) ||
        throw(ArgumentError("$(label) diagnostic buffer must match the wavefront size"))
    eltype(buffer) === T ||
        throw(ArgumentError("$(label) diagnostic buffer must use wavefront precision $(T)"))
    Proper.same_backend_style(typeof(buffer), typeof(wfo.field)) ||
        throw(ArgumentError("$(label) diagnostic buffer must use the wavefront array backend"))
    return nothing
end

function _capture_after_occulter!(diagnostics::CoronagraphDiagnostics, wfo::WaveFront)
    _check_coronagraph_diagnostic_buffer(
        diagnostics.after_occulter_amplitude,
        wfo,
        "after-occulter",
    )
    copyto!(diagnostics.after_occulter_amplitude, prop_get_amplitude(wfo))
    return nothing
end

function _capture_before_lyot!(diagnostics::CoronagraphDiagnostics, wfo::WaveFront)
    _check_coronagraph_diagnostic_buffer(
        diagnostics.before_lyot_amplitude,
        wfo,
        "before-Lyot",
    )
    copyto!(diagnostics.before_lyot_amplitude, prop_get_amplitude(wfo))
    return nothing
end

function _apply_occulter!(wfo::WaveFront, ::GaussianOcculter, occrad::Real, occrad_m::Real)
    r = prop_radius(wfo)
    h = sqrt(-0.5 * occrad_m^2 / log(1 - sqrt(0.5)))
    gauss_spot = @. 1 - exp(-0.5 * (r / h)^2)
    prop_multiply(wfo, gauss_spot)
    return wfo
end

function _apply_occulter!(wfo::WaveFront, ::SolidOcculter, occrad::Real, occrad_m::Real)
    prop_circular_obscuration(wfo, occrad_m)
    return wfo
end

function _apply_occulter!(wfo::WaveFront, ::EighthOrderOcculter, occrad::Real, occrad_m::Real)
    # `prop_8th_order_mask` interprets `occrad` in lambda/D units here.
    prop_8th_order_mask(wfo, occrad; circular=true)
    return wfo
end

function _apply_lyot_stop!(wfo::WaveFront, ::GaussianOcculter)
    prop_circular_aperture(wfo, 0.25; NORM=true)
    return wfo
end

function _apply_lyot_stop!(wfo::WaveFront, ::SolidOcculter)
    prop_circular_aperture(wfo, 0.84; NORM=true)
    return wfo
end

function _apply_lyot_stop!(wfo::WaveFront, ::EighthOrderOcculter)
    prop_circular_aperture(wfo, 0.50; NORM=true)
    return wfo
end

function coronagraph(
    wfo::WaveFront,
    f_lens::Real,
    occulter,
    diam::Real;
    diagnostics::Union{Nothing,CoronagraphDiagnostics}=nothing,
)
    occulter = normalize_occulter(occulter)

    prop_lens(wfo, f_lens, "coronagraph imaging lens")
    prop_propagate(wfo, f_lens, "occulter")

    λ = prop_get_wavelength(wfo)
    occrad = 4.0
    occrad_rad = occrad * λ / diam
    dx_m = prop_get_sampling(wfo)
    dx_rad = prop_get_sampling_radians(wfo)
    occrad_m = occrad_rad * dx_m / dx_rad

    _apply_occulter!(wfo, occulter, occrad, occrad_m)
    _capture_after_occulter!(diagnostics, wfo)

    prop_propagate(wfo, f_lens, "pupil reimaging lens")
    prop_lens(wfo, f_lens, "pupil reimaging lens")
    prop_propagate(wfo, 2f_lens, "lyot stop")
    _capture_before_lyot!(diagnostics, wfo)

    _apply_lyot_stop!(wfo, occulter)

    prop_propagate(wfo, f_lens, "reimaging lens")
    prop_lens(wfo, f_lens, "reimaging lens")
    prop_propagate(wfo, f_lens, "final focus")
    return wfo
end
