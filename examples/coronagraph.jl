using Proper
using Plots

function coronagraph(wfo::WaveFront, f_lens::Real, occulter_type::AbstractString, diam::Real)
    prop_lens(wfo, f_lens, "coronagraph imaging lens")
    prop_propagate(wfo, f_lens, "occulter")

    λ = prop_get_wavelength(wfo)
    occrad = 4.0
    occrad_rad = occrad * λ / diam
    dx_m = prop_get_sampling(wfo)
    dx_rad = prop_get_sampling_radians(wfo)
    occrad_m = occrad_rad * dx_m / dx_rad

    if occulter_type == "GAUSSIAN"
        r = prop_radius(wfo)
        h = sqrt(-0.5 * occrad_m^2 / log(1 - sqrt(0.5)))
        gauss_spot = @. 1 - exp(-0.5 * (r / h)^2)
        prop_multiply(wfo, gauss_spot)
    elseif occulter_type == "SOLID"
        prop_circular_obscuration(wfo, occrad_m)
    elseif occulter_type == "8TH_ORDER"
        prop_8th_order_mask(wfo, occrad; circular=true)
    else
        throw(ArgumentError("Unknown occulter_type=$occulter_type"))
    end

    prop_propagate(wfo, f_lens, "pupil reimaging lens")
    prop_lens(wfo, f_lens, "pupil reimaging lens")
    prop_propagate(wfo, 2f_lens, "lyot stop")

    if occulter_type == "GAUSSIAN"
        prop_circular_aperture(wfo, 0.25; NORM=true)
    elseif occulter_type == "SOLID"
        prop_circular_aperture(wfo, 0.84; NORM=true)
    else
        prop_circular_aperture(wfo, 0.50; NORM=true)
    end

    prop_propagate(wfo, f_lens, "reimaging lens")
    prop_lens(wfo, f_lens, "reimaging lens")
    prop_propagate(wfo, f_lens, "final focus")
    return wfo
end
