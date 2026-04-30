using Proper
using Plots
include(joinpath(@__DIR__, "_shared.jl"))
include(joinpath(@__DIR__, "_passvalue.jl"))

function hubble_simple(wavelength::Real, gridsize::Integer, passvalue; kwargs...)
    return hubble_simple(wavelength, gridsize; passvalue_kwargs(passvalue)..., kwargs...)
end

function hubble_simple(wavelength::Real, gridsize::Integer; delta_sec::Real=0.0)
    diam = 2.4
    fl_pri = 5.52085
    d_pri_sec = 4.907028205
    fl_sec = -0.6790325
    d_sec_to_focus = 6.3919974
    beam_ratio = 0.5

    wfo = prop_begin(diam, wavelength, gridsize; beam_diam_fraction=beam_ratio)
    prop_circular_aperture(wfo, diam / 2)
    prop_circular_obscuration(wfo, 0.396)
    prop_rectangular_obscuration(wfo, 0.0264, 2.5)
    prop_rectangular_obscuration(wfo, 2.5, 0.0264)
    prop_circular_obscuration(wfo, 0.078, -0.9066, -0.5538)
    prop_circular_obscuration(wfo, 0.078, 0.0, 1.0705)
    prop_circular_obscuration(wfo, 0.078, 0.9127, -0.5477)

    prop_define_entrance(wfo)
    prop_lens(wfo, fl_pri, "primary")
    prop_propagate(wfo, d_pri_sec + delta_sec, "secondary")
    prop_lens(wfo, fl_sec, "secondary")
    prop_propagate(wfo, d_sec_to_focus + delta_sec, "HST focus"; TO_PLANE=false)
    return prop_end(wfo)
end

if abspath(PROGRAM_FILE) == @__FILE__
    psf, sampling = hubble_simple(0.55e-6, 512)
    println("hubble_simple: sampling = ", sampling)
    plot_psf(psf; title="hubble_simple")
end
