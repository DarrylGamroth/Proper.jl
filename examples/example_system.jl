using proper
using Plots
include(joinpath(@__DIR__, "_shared.jl"))

function example_system(wavelength::Real, gridsize::Integer)
    diam = 1.0
    lens_fl = 20.0
    beam_ratio = 0.5

    wfo = prop_begin(diam, wavelength, gridsize; beam_diam_fraction=beam_ratio)

    if !prop_is_statesaved(wfo)
        prop_circular_aperture(wfo, diam / 2)
        prop_define_entrance(wfo)
        prop_lens(wfo, lens_fl, "1st lens")
        prop_propagate(wfo, lens_fl, "intermediate focus")
    end

    prop_state(wfo)
    prop_propagate(wfo, lens_fl, "second lens")
    prop_lens(wfo, lens_fl, "second lens")
    return prop_end(wfo)
end

if abspath(PROGRAM_FILE) == @__FILE__
    psf, sampling = example_system(0.55e-6, 256)
    println("example_system: sampling = ", sampling)
    plot_psf(psf; title="example_system")
end
