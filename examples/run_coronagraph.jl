using proper
include(joinpath(@__DIR__, "telescope.jl"))
include(joinpath(@__DIR__, "coronagraph.jl"))

function run_coronagraph(wavelength::Real, grid_size::Integer, passvalue=Dict("use_errors" => false, "occulter_type" => "GAUSSIAN"))
    diam = 0.1
    f_lens = 24 * diam
    beam_ratio = 0.3

    wfo = prop_begin(diam, wavelength, grid_size; beam_diam_fraction=beam_ratio)
    prop_circular_aperture(wfo, diam / 2)
    prop_define_entrance(wfo)

    telescope(wfo, f_lens, get(passvalue, "use_errors", get(passvalue, :use_errors, false)))
    coronagraph(wfo, f_lens, get(passvalue, "occulter_type", get(passvalue, :occulter_type, "GAUSSIAN")), diam)

    return prop_end(wfo)
end
