Base.@kwdef struct ProperConfig{P<:CompatPolicy}
    policy::P = Python334Policy()
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
end

function WaveFront(field::A, wavelength_m::T, sampling_m::T, z_m::T, beam_diameter_m::T) where {T<:AbstractFloat,A<:AbstractMatrix{Complex{T}}}
    return WaveFront{T,A}(field, wavelength_m, sampling_m, z_m, beam_diameter_m)
end
