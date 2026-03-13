Base.@kwdef struct ProperConfig
    verbose::Bool = false
end

mutable struct ProperRuntime
    run_id::Int
end

ProperRuntime() = ProperRuntime(0)

@enum ReferenceSurface::UInt8 PLANAR SPHERICAL
@enum BeamType::UInt8 INSIDE OUTSIDE
@enum PropagatorType::UInt8 INSIDE_TO_INSIDE INSIDE_TO_OUTSIDE OUTSIDE_TO_INSIDE OUTSIDE_TO_OUTSIDE

@inline function propagator_transition(old::BeamType, new::BeamType)::PropagatorType
    if old === INSIDE
        return new === INSIDE ? INSIDE_TO_INSIDE : INSIDE_TO_OUTSIDE
    end
    return new === INSIDE ? OUTSIDE_TO_INSIDE : OUTSIDE_TO_OUTSIDE
end

@inline function to_plane_propagator(pt::PropagatorType)::PropagatorType
    if pt === INSIDE_TO_OUTSIDE
        return INSIDE_TO_INSIDE
    elseif pt === OUTSIDE_TO_OUTSIDE
        return OUTSIDE_TO_INSIDE
    end
    return pt
end

mutable struct WaveFront{T,A<:AbstractMatrix{Complex{T}},WS<:ProperWorkspace{T}}
    field::A
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
    workspace::WS
end

@inline function field_backend_template(ws::ProperWorkspace)
    return ws.fft.scratch
end

@inline workspace_float_type(::Nothing) = nothing
@inline workspace_float_type(ws::ProperWorkspace{T}) where {T<:AbstractFloat} = T

@inline function unit_field(::Type{T}, n::Integer, ::Nothing) where {T<:AbstractFloat}
    return fill(complex(one(T), zero(T)), n, n)
end

@inline function unit_field(::Type{T}, n::Integer, ws::ProperWorkspace) where {T<:AbstractFloat}
    field = similar(field_backend_template(ws), Complex{T}, n, n)
    fill!(field, complex(one(T), zero(T)))
    return field
end

function WaveFront(
    field::A,
    wavelength_m::T,
    sampling_m::T,
    z_m::T,
    beam_diameter_m::T,
    ws::WS,
) where {T<:AbstractFloat,A<:AbstractMatrix{Complex{T}},WS<:ProperWorkspace{T}}
    w0_m = beam_diameter_m / 2
    z_ray = pi * w0_m^2 / wavelength_m
    return WaveFront(
        field,
        wavelength_m,
        sampling_m,
        z_m,
        beam_diameter_m,
        zero(T),
        w0_m,
        z_ray,
        T(1e9),
        PLANAR,
        INSIDE,
        INSIDE_TO_INSIDE,
        one(T),
        ws,
    )
end

function WaveFront(field::A, wavelength_m::T, sampling_m::T, z_m::T, beam_diameter_m::T) where {T<:AbstractFloat,A<:AbstractMatrix{Complex{T}}}
    return WaveFront(field, wavelength_m, sampling_m, z_m, beam_diameter_m, ProperWorkspace(A, T))
end
