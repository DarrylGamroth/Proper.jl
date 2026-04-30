using Proper
using Plots
include(joinpath(@__DIR__, "_passvalue.jl"))

function multi_example(lambda_m::Real, n::Integer, passvalue; kwargs...)
    return multi_example(lambda_m, n; passvalue_kwargs(passvalue)..., kwargs...)
end

function multi_example(lambda_m::Real, n::Integer; use_dm::Bool=false, dm=nothing)
    diam = 0.048
    pupil_ratio = 0.25
    fl_lens = 0.48

    wfo = prop_begin(diam, lambda_m, n; beam_diam_fraction=pupil_ratio)
    prop_circular_aperture(wfo, diam / 2)
    prop_define_entrance(wfo)

    if use_dm
        prop_dm(wfo, dm === nothing ? zeros(n, n) : dm; mirror=false)
    end

    prop_lens(wfo, fl_lens)
    prop_propagate(wfo, fl_lens)
    return prop_end(wfo; noabs=true)
end

if abspath(PROGRAM_FILE) == @__FILE__
    model = prepare_model(
        multi_example,
        0.55,
        256;
        name=:multi_example,
        PASSVALUE=[(; use_dm=false), (; use_dm=true, dm=zeros(256, 256))],
        pool_size=2,
    )
    stack, samplings = prop_run_multi(model)
    println("multi_example: samplings = ", samplings)
    heatmap(log10.(abs.(stack[:, :, 1]) .+ eps()); aspect_ratio=:equal, title="multi_example [0]")
end
