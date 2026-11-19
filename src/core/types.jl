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

"""
    WaveFrontGeometry(diameter_m, wavelength_m, gridsize; beam_diam_fraction=0.5)

Prepare the scalar geometry of a square PROPER wavefront at its entrance pupil.
The floating-point wavelength type determines the geometry precision. Advance
the geometry with `prop_lens` and `prop_propagate`. The grid size participates
in sampling and aperture validation. This object owns no electric field or
FFT workspace and is intended for preparing a fixed prescription.
"""
mutable struct WaveFrontGeometry{T<:AbstractFloat}
    const gridsize::Int
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

    function WaveFrontGeometry(diameter_m::Real, wavelength_m::Real,
            gridsize::Integer; beam_diam_fraction::Real=0.5)
        n = Int(gridsize)
        n > 0 || throw(ArgumentError("gridsize must be positive"))
        T = float(typeof(wavelength_m))
        λ, d = T(wavelength_m), T(diameter_m)
        sampling = _entrance_sampling(T, d, n, beam_diam_fraction)
        return new{T}(n, _entrance_geometry(λ, sampling, zero(T), d)...)
    end
end

const OpticalSurfaceState = Union{WaveFront,WaveFrontGeometry}

@inline _entrance_sampling(::Type{T}, diameter, n, fraction) where {T} =
    diameter / (T(n) * T(fraction))

@inline function _entrance_geometry(wavelength::T, sampling::T, z::T, diameter::T) where {T}
    waist = diameter / 2
    rayleigh = pi * waist^2 / wavelength
    return (wavelength, sampling, z, diameter, zero(T), waist, rayleigh,
        T(1e9), PLANAR, INSIDE, INSIDE_TO_INSIDE, one(T))
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
    return WaveFront(
        field,
        _entrance_geometry(wavelength_m, sampling_m, z_m, beam_diameter_m)...,
        ws,
    )
end

@inline _require_square_propagation_grid(wf::WaveFrontGeometry, ::Symbol) = wf.gridsize

function WaveFront(field::A, wavelength_m::T, sampling_m::T, z_m::T, beam_diameter_m::T) where {T<:AbstractFloat,A<:AbstractMatrix{Complex{T}}}
    return WaveFront(field, wavelength_m, sampling_m, z_m, beam_diameter_m, ProperWorkspace(A, T))
end

@inline function _require_square_propagation_grid(wf::WaveFront, routine::Symbol)
    ny, nx = size(wf.field)
    ny > 0 && ny == nx || throw(ArgumentError(
        "$(routine) requires a non-empty square wavefront field; got $(ny)x$(nx)",
    ))
    return ny
end
