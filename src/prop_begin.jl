"""Initialize a wavefront at entrance pupil."""
function prop_begin(diam::Real, wavelength_m::Real, gridsize::Integer; beam_diam_fraction::Real=0.5)
    n = Int(gridsize)
    n > 0 || throw(ArgumentError("gridsize must be positive"))
    λ = float(wavelength_m)
    d = float(diam)
    ndiam = n * float(beam_diam_fraction)
    sampling = d / ndiam
    field = fill(complex(one(λ), zero(λ)), n, n)
    return WaveFront(field, λ, sampling, zero(λ), d)
end
