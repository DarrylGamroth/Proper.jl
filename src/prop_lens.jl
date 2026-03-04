"""Apply a thin lens with focal length `lens_fl` in meters."""
function prop_lens(wf::WaveFront, lens_fl::Real, surface_name::AbstractString="")
    iszero(lens_fl) && throw(ArgumentError("lens_fl must be non-zero"))
    r = prop_radius(wf)
    phase = @. -pi * (r^2) / (wf.wavelength_m * lens_fl)
    wf.field .*= cis.(phase)
    wf.z_m = wf.z_m
    return wf
end
