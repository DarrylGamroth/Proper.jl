"""Initialize a wavefront at entrance pupil."""
function prop_begin(diam::Real, wavelength_m::Real, gridsize::Integer; beam_diam_fraction::Real=1.0)
    n = Int(gridsize)
    n > 0 || throw(ArgumentError("gridsize must be positive"))
    d = float(diam) * float(beam_diam_fraction)
    return prop_wavefront(n, float(wavelength_m), d)
end
