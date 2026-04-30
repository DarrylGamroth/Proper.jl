using Proper
using Plots

function _coronagraph_plot_title(occulter_type::AbstractString)
    if occulter_type == "GAUSSIAN"
        return "Gaussian spot"
    elseif occulter_type == "SOLID"
        return "Solid spot"
    elseif occulter_type == "8TH_ORDER"
        return "8th order band limited spot"
    else
        return "Coronagraph"
    end
end

function _plot_coronagraph_planes(after_occulter::AbstractMatrix, before_lyot::AbstractMatrix, occulter_type::AbstractString)
    p1 = heatmap(
        after_occulter;
        aspect_ratio=:equal,
        color=:grays,
        colorbar=true,
        title="After Occulter",
    )
    p2 = heatmap(
        before_lyot;
        aspect_ratio=:equal,
        color=:grays,
        colorbar=true,
        title="Before Lyot Stop",
    )
    plt = plot(p1, p2; layout=(1, 2), size=(1200, 500), plot_title=_coronagraph_plot_title(occulter_type))
    display(plt)
    return plt
end

function coronagraph(wfo::WaveFront, f_lens::Real, occulter_type::AbstractString, diam::Real; PLOT::Bool=true)
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
        # `prop_8th_order_mask` interprets `occrad` in lambda/D units here.
        prop_8th_order_mask(wfo, occrad; circular=true)
    else
        throw(ArgumentError("Unknown occulter_type=$occulter_type"))
    end

    after_occulter = PLOT ? sqrt.(prop_get_amplitude(wfo)) : nothing

    prop_propagate(wfo, f_lens, "pupil reimaging lens")
    prop_lens(wfo, f_lens, "pupil reimaging lens")
    prop_propagate(wfo, 2f_lens, "lyot stop")

    if PLOT
        before_lyot = prop_get_amplitude(wfo) .^ 0.2
        _plot_coronagraph_planes(after_occulter, before_lyot, occulter_type)
    end

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
