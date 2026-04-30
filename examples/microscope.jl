using Proper
using Plots
include(joinpath(@__DIR__, "_shared.jl"))
include(joinpath(@__DIR__, "_passvalue.jl"))

function microscope(wavelength::Real, gridsize::Integer, passvalue; kwargs...)
    return microscope(wavelength, gridsize; passvalue_kwargs(passvalue)..., kwargs...)
end

function microscope(wavelength::Real, gridsize::Integer; focus_offset::Real=0.0)
    d_objective = 0.005
    fl_objective = 0.010
    fl_eyepiece = 0.020
    fl_eye = 0.022
    beam_ratio = 0.4

    wfo = prop_begin(d_objective, wavelength, gridsize; beam_diam_fraction=beam_ratio)

    d1 = 0.160
    d_intermediate_image = fl_objective + d1
    d_object = 1 / (1 / fl_objective - 1 / d_intermediate_image)

    prop_circular_aperture(wfo, d_objective / 2)
    prop_define_entrance(wfo)

    prop_lens(wfo, -(d_object + focus_offset))
    prop_lens(wfo, fl_objective, "objective")

    prop_propagate(wfo, d_intermediate_image, "intermediate image")
    prop_propagate(wfo, fl_eyepiece, "eyepiece")
    prop_lens(wfo, fl_eyepiece, "eyepiece")
    exit_pupil_distance = fl_eyepiece / (1 - fl_eyepiece / (d_intermediate_image + fl_eyepiece))
    prop_propagate(wfo, exit_pupil_distance, "exit pupil/eye")
    prop_lens(wfo, fl_eye, "eye")
    prop_propagate(wfo, fl_eye, "retina")

    return prop_end(wfo)
end

if abspath(PROGRAM_FILE) == @__FILE__
    psf, sampling = microscope(0.55e-6, 256)
    println("microscope: sampling = ", sampling)
    plot_psf(psf; title="microscope")
end
