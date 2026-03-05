using Proper
using Plots
include(joinpath(@__DIR__, "_shared.jl"))

function simple_telescope(wavelength::Real, gridsize::Integer)
    d_objective = 0.060
    fl_objective = 15.0 * d_objective
    fl_eyepiece = 0.021
    fl_eye = 0.022
    beam_ratio = 0.5

    wfo = prop_begin(d_objective, wavelength, gridsize; beam_diam_fraction=beam_ratio)
    prop_circular_aperture(wfo, d_objective / 2)
    prop_define_entrance(wfo)
    prop_lens(wfo, fl_objective, "objective")
    prop_propagate(wfo, fl_objective + fl_eyepiece, "eyepiece")
    prop_lens(wfo, fl_eyepiece, "eyepiece")

    exit_pupil_distance = fl_eyepiece / (1 - fl_eyepiece / (fl_objective + fl_eyepiece))
    prop_propagate(wfo, exit_pupil_distance, "exit pupil at eye lens")
    prop_lens(wfo, fl_eye, "eye")
    prop_propagate(wfo, fl_eye, "retina")

    return prop_end(wfo)
end

if abspath(PROGRAM_FILE) == @__FILE__
    psf, sampling = simple_telescope(0.55e-6, 256)
    println("simple_telescope: sampling = ", sampling)
    plot_psf(psf; title="simple_telescope")
end
