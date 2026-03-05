"""Create a wavefront container with a complex field initialized to 1+0im."""
function prop_wavefront(gridsize::Integer, wavelength_m::Real, beam_diameter_m::Real; sampling_m::Union{Nothing,Real}=nothing)
    n = Int(gridsize)
    n > 0 || throw(ArgumentError("gridsize must be positive"))
    λ = float(wavelength_m)
    d = float(beam_diameter_m)
    s = sampling_m === nothing ? d / n : float(sampling_m)
    field = fill(complex(one(λ), zero(λ)), n, n)
    w0 = d / 2
    zray = pi * w0^2 / λ
    return WaveFront(
        field,
        λ,
        s,
        zero(λ),
        d,
        zero(λ),
        w0,
        zray,
        float(1e9),
        PLANAR,
        INSIDE_,
        INSIDE__to_INSIDE_,
        one(λ),
    )
end
