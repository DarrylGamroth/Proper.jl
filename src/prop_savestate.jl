using Serialization

struct WaveFrontState{T<:AbstractFloat}
    field::Matrix{Complex{T}}
    wavelength_m::T
    sampling_m::T
    z_m::T
    beam_diameter_m::T
    z_w0_m::T
    w0_m::T
    z_rayleigh_m::T
    current_fratio::T
    reference_surface::ReferenceSurface
    beam_type_old::BeamType
    propagator_type::PropagatorType
    rayleigh_factor::T
end

@inline function WaveFrontState(wf::WaveFront{T}) where {T<:AbstractFloat}
    return WaveFrontState{T}(
        Matrix{Complex{T}}(wf.field),
        wf.wavelength_m,
        wf.sampling_m,
        wf.z_m,
        wf.beam_diameter_m,
        wf.z_w0_m,
        wf.w0_m,
        wf.z_rayleigh_m,
        wf.current_fratio,
        wf.reference_surface,
        wf.beam_type_old,
        wf.propagator_type,
        wf.rayleigh_factor,
    )
end

"""Serialize wavefront state to disk."""
function prop_savestate(wf::WaveFront, path::AbstractString)
    open(path, "w") do io
        serialize(io, WaveFrontState(wf))
    end
    return path
end
