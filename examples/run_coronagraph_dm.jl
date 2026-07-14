using Proper
using Random
include(joinpath(@__DIR__, "_passvalue.jl"))
include(joinpath(@__DIR__, "telescope_dm.jl"))
include(joinpath(@__DIR__, "coronagraph.jl"))

function run_coronagraph_dm(wavelength::Real, grid_size::Integer, passvalue; kwargs...)
    return run_coronagraph_dm(wavelength, grid_size; passvalue_kwargs(passvalue)..., kwargs...)
end

function run_coronagraph_dm(
    wavelength::Real,
    grid_size::Integer;
    use_errors::Bool=false,
    use_dm::Bool=false,
    occulter=:gaussian,
    occulter_type=nothing,
    diagnostics::Union{Nothing,CoronagraphDiagnostics}=nothing,
    map_file::AbstractString="telescope_obj.fits",
    rng::AbstractRNG=Random.default_rng(),
)
    diam = 0.1
    f_lens = 24 * diam
    beam_ratio = 0.3

    wfo = prop_begin(diam, wavelength, grid_size; beam_diam_fraction=beam_ratio)
    prop_circular_aperture(wfo, diam / 2)
    prop_define_entrance(wfo)

    telescope_dm(
        wfo,
        f_lens,
        use_errors,
        use_dm,
        ;
        map_file,
        rng,
    )
    coronagraph(
        wfo,
        f_lens,
        occulter_type === nothing ? occulter : occulter_type,
        diam;
        diagnostics,
    )

    return prop_end(wfo)
end
