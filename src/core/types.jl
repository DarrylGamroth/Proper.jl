Base.@kwdef struct ProperConfig
    verbose::Bool = false
end

mutable struct ProperRuntime
    run_id::Int
end

ProperRuntime() = ProperRuntime(0)

mutable struct WaveFront{T,A<:AbstractMatrix{Complex{T}}}
    field::A
    wavelength_m::T
    sampling_m::T
    z_m::T
    beam_diameter_m::T
    z_w0_m::T
    w0_m::T
    z_rayleigh_m::T
    current_fratio::T
    reference_surface::Symbol
    beam_type_old::Symbol
    propagator_type::Symbol
    rayleigh_factor::T
end

function WaveFront(
    field::A,
    wavelength_m::T,
    sampling_m::T,
    z_m::T,
    beam_diameter_m::T,
    z_w0_m::T,
    w0_m::T,
    z_rayleigh_m::T,
    current_fratio::T,
    reference_surface::Symbol,
    beam_type_old::Symbol,
    propagator_type::Symbol,
    rayleigh_factor::T,
) where {T<:AbstractFloat,A<:AbstractMatrix{Complex{T}}}
    return WaveFront{T,A}(
        field,
        wavelength_m,
        sampling_m,
        z_m,
        beam_diameter_m,
        z_w0_m,
        w0_m,
        z_rayleigh_m,
        current_fratio,
        reference_surface,
        beam_type_old,
        propagator_type,
        rayleigh_factor,
    )
end

function WaveFront(field::A, wavelength_m::T, sampling_m::T, z_m::T, beam_diameter_m::T) where {T<:AbstractFloat,A<:AbstractMatrix{Complex{T}}}
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
        :PLANAR,
        :INSIDE_,
        :INSIDE__to_INSIDE_,
        one(T),
    )
end
