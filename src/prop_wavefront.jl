"""Create a wavefront container with a complex field initialized to 1+0im."""
function prop_wavefront(gridsize::Integer, wavelength_m::Real, beam_diameter_m::Real; sampling_m::Union{Nothing,Real}=nothing)
    n = Int(gridsize)
    n > 0 || throw(ArgumentError("gridsize must be positive"))
    λ = float(wavelength_m)
    d = float(beam_diameter_m)
    s = sampling_m === nothing ? d / n : float(sampling_m)
    field = fill(complex(one(λ), zero(λ)), n, n)
    return WaveFront(field, λ, s, zero(λ), d)
end
