using Proper
include(joinpath(@__DIR__, "_passvalue.jl"))

function example_system(wavelength::Real, gridsize::Integer, passvalue; kwargs...)
    return example_system(wavelength, gridsize; passvalue_kwargs(passvalue)..., kwargs...)
end

function example_system(wavelength::Real, gridsize::Integer; state_path=tempname())
    diam = 1.0
    lens_fl = 20.0
    beam_ratio = 0.5

    wfo = prop_begin(diam, wavelength, gridsize; beam_diam_fraction=beam_ratio)

    if !prop_is_statesaved(state_path)
        prop_circular_aperture(wfo, diam / 2)
        prop_define_entrance(wfo)
        prop_lens(wfo, lens_fl, "1st lens")
        prop_propagate(wfo, lens_fl, "intermediate focus")
        prop_savestate(wfo, state_path)
    else
        prop_state(wfo, state_path)
    end

    prop_propagate(wfo, lens_fl, "second lens")
    prop_lens(wfo, lens_fl, "second lens")
    return prop_end(wfo)
end

if abspath(PROGRAM_FILE) == @__FILE__
    using Plots

    psf, sampling = example_system(0.55e-6, 256)
    println("example_system: sampling = ", sampling)
    display(heatmap(psf .^ 0.25; aspect_ratio=:equal, color=:grays, title="example_system"))
end
