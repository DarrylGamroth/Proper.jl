using Proper
using Plots
include(joinpath(@__DIR__, "_shared.jl"))
include(joinpath(@__DIR__, "_passvalue.jl"))

function psdtest(wavelength::Real, gridsize::Integer, passvalue; kwargs...)
    return psdtest(wavelength, gridsize; passvalue_kwargs(passvalue)..., kwargs...)
end

function psdtest(wavelength::Real, gridsize::Integer; usepsdmap::Bool=true)
    lens_diam = 0.212
    lens_fl = 24 * lens_diam
    beam_width_ratio = 0.5

    wfo = prop_begin(lens_diam, wavelength, gridsize; beam_diam_fraction=beam_width_ratio)
    prop_circular_aperture(wfo, lens_diam / 2)
    prop_define_entrance(wfo)

    if usepsdmap
        a = 3.29e-23
        b = 212.26
        c = 7.8
        prop_psd_errormap(wfo, a, b, c)
    else
        # This is a wavefront error map in nanometers, so apply it explicitly as a
        # wavefront map after converting nm -> m.
        prop_errormap(wfo, "errormap.fits"; SAMPLING=0.0004, MULTIPLY=1e-9, WAVEFRONT=true)
    end

    prop_lens(wfo, lens_fl, "telescope lens")
    prop_propagate(wfo, prop_get_distancetofocus(wfo), "intermediate focus")
    # HWHM is in lambda/D units unless `METERS=true` is supplied.
    prop_8th_order_mask(wfo, 4; circular=true)

    prop_propagate(wfo, lens_fl, "pupil imaging lens")
    prop_lens(wfo, lens_fl, "pupil imaging lens")
    prop_propagate(wfo, lens_fl, "lyot stop")
    prop_circular_aperture(wfo, 0.53; NORM=true)

    prop_propagate(wfo, lens_fl, "reimaging lens")
    prop_lens(wfo, lens_fl, "reimaging lens")
    prop_propagate(wfo, prop_get_distancetofocus(wfo), "final focus")

    return prop_end(wfo)
end

if abspath(PROGRAM_FILE) == @__FILE__
    psf, sampling = psdtest(0.55e-6, 256)
    println("psdtest: sampling = ", sampling)
    plot_psf(psf; title="psdtest")
end
