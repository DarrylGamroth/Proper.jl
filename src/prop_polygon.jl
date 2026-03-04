"""Regular polygon aperture approximation using circumscribed radius."""
function prop_polygon(wf::WaveFront, nsides::Integer, radius::Real, xc::Real=0.0, yc::Real=0.0; kwargs...)
    nsides >= 3 || throw(ArgumentError("nsides must be >= 3"))
    # Placeholder: circle approximation keeps API stable until exact polygon kernel lands.
    return prop_ellipse(wf, radius, radius, xc, yc; kwargs...)
end
