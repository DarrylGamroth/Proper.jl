Base.@kwdef struct ProperConfig
    verbose::Bool = false
end

mutable struct ProperRuntime
    run_id::Int
end

ProperRuntime() = ProperRuntime(0)

@enum ReferenceSurface::UInt8 PLANAR SPHERI
@enum BeamType::UInt8 INSIDE_ OUTSIDE
@enum PropagatorType::UInt8 INSIDE__to_INSIDE_ INSIDE__to_OUTSIDE OUTSIDE_to_INSIDE_ OUTSIDE_to_OUTSIDE

@inline reference_surface(v::ReferenceSurface) = v
@inline beam_type(v::BeamType) = v
@inline propagator_type(v::PropagatorType) = v

@inline function reference_surface(v::Symbol)
    v === :PLANAR && return PLANAR
    v === :SPHERI && return SPHERI
    throw(ArgumentError("Unknown reference surface symbol: $v"))
end

@inline function beam_type(v::Symbol)
    v === :INSIDE_ && return INSIDE_
    v === :OUTSIDE && return OUTSIDE
    throw(ArgumentError("Unknown beam type symbol: $v"))
end

@inline function propagator_type(v::Symbol)
    v === :INSIDE__to_INSIDE_ && return INSIDE__to_INSIDE_
    v === :INSIDE__to_OUTSIDE && return INSIDE__to_OUTSIDE
    v === :OUTSIDE_to_INSIDE_ && return OUTSIDE_to_INSIDE_
    v === :OUTSIDE_to_OUTSIDE && return OUTSIDE_to_OUTSIDE
    throw(ArgumentError("Unknown propagator type symbol: $v"))
end

@inline reference_surface(v::AbstractString) = reference_surface(Symbol(v))
@inline beam_type(v::AbstractString) = beam_type(Symbol(v))
@inline propagator_type(v::AbstractString) = propagator_type(Symbol(v))

@inline reference_surface_symbol(v::ReferenceSurface)::Symbol = (v === PLANAR ? :PLANAR : :SPHERI)
@inline beam_type_symbol(v::BeamType)::Symbol = (v === INSIDE_ ? :INSIDE_ : :OUTSIDE)

@inline function propagator_type_symbol(v::PropagatorType)::Symbol
    v === INSIDE__to_INSIDE_ && return :INSIDE__to_INSIDE_
    v === INSIDE__to_OUTSIDE && return :INSIDE__to_OUTSIDE
    v === OUTSIDE_to_INSIDE_ && return :OUTSIDE_to_INSIDE_
    return :OUTSIDE_to_OUTSIDE
end

Base.convert(::Type{ReferenceSurface}, v::ReferenceSurface) = v
Base.convert(::Type{BeamType}, v::BeamType) = v
Base.convert(::Type{PropagatorType}, v::PropagatorType) = v

Base.convert(::Type{ReferenceSurface}, v::Symbol) = reference_surface(v)
Base.convert(::Type{BeamType}, v::Symbol) = beam_type(v)
Base.convert(::Type{PropagatorType}, v::Symbol) = propagator_type(v)

Base.convert(::Type{ReferenceSurface}, v::AbstractString) = reference_surface(v)
Base.convert(::Type{BeamType}, v::AbstractString) = beam_type(v)
Base.convert(::Type{PropagatorType}, v::AbstractString) = propagator_type(v)

@inline Base.:(==)(lhs::ReferenceSurface, rhs::Symbol) = reference_surface_symbol(lhs) === rhs
@inline Base.:(==)(lhs::BeamType, rhs::Symbol) = beam_type_symbol(lhs) === rhs
@inline Base.:(==)(lhs::PropagatorType, rhs::Symbol) = propagator_type_symbol(lhs) === rhs
@inline Base.:(==)(lhs::Symbol, rhs::ReferenceSurface) = rhs == lhs
@inline Base.:(==)(lhs::Symbol, rhs::BeamType) = rhs == lhs
@inline Base.:(==)(lhs::Symbol, rhs::PropagatorType) = rhs == lhs

@inline function propagator_transition(old::BeamType, new::BeamType)::PropagatorType
    if old === INSIDE_
        return new === INSIDE_ ? INSIDE__to_INSIDE_ : INSIDE__to_OUTSIDE
    end
    return new === INSIDE_ ? OUTSIDE_to_INSIDE_ : OUTSIDE_to_OUTSIDE
end

@inline function to_plane_propagator(pt::PropagatorType)::PropagatorType
    if pt === INSIDE__to_OUTSIDE
        return INSIDE__to_INSIDE_
    elseif pt === OUTSIDE_to_OUTSIDE
        return OUTSIDE_to_INSIDE_
    end
    return pt
end

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
    reference_surface::ReferenceSurface
    beam_type_old::BeamType
    propagator_type::PropagatorType
    rayleigh_factor::T
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
        PLANAR,
        INSIDE_,
        INSIDE__to_INSIDE_,
        one(T),
    )
end
