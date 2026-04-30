using Proper
include(joinpath(@__DIR__, "_passvalue.jl"))
include(joinpath(@__DIR__, "telescope.jl"))
include(joinpath(@__DIR__, "coronagraph.jl"))

function run_coronagraph(wavelength::Real, grid_size::Integer, passvalue; kwargs...)
    return run_coronagraph(wavelength, grid_size; passvalue_kwargs(passvalue)..., kwargs...)
end

function run_coronagraph(
    wavelength::Real,
    grid_size::Integer;
    use_errors::Bool=false,
    occulter=:gaussian,
    occulter_type=nothing,
    plot::Bool=false,
)
    diam = 0.1
    f_lens = 24 * diam
    beam_ratio = 0.3

    wfo = prop_begin(diam, wavelength, grid_size; beam_diam_fraction=beam_ratio)
    prop_circular_aperture(wfo, diam / 2)
    prop_define_entrance(wfo)

    telescope(wfo, f_lens, use_errors)
    coronagraph(
        wfo,
        f_lens,
        occulter_type === nothing ? occulter : occulter_type,
        diam;
        plot=plot,
    )

    return prop_end(wfo)
end
