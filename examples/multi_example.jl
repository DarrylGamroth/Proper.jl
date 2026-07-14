using Proper
include(joinpath(@__DIR__, "_passvalue.jl"))

function multi_example(lambda_m::Real, n::Integer, passvalue; kwargs...)
    return multi_example(lambda_m, n; passvalue_kwargs(passvalue)..., kwargs...)
end

function multi_example(lambda_m::Real, n::Integer; use_dm::Bool=false, dm=nothing)
    diam = 0.048
    pupil_ratio = 0.25
    fl_lens = 0.48
    n_actuators = 48

    wfo = prop_begin(diam, lambda_m, n; beam_diam_fraction=pupil_ratio)
    prop_circular_aperture(wfo, diam / 2)
    prop_define_entrance(wfo)

    if use_dm
        dm_surface = dm === nothing ? zeros(typeof(wfo.sampling_m), n_actuators, n_actuators) : dm
        size(dm_surface) == (n_actuators, n_actuators) || throw(ArgumentError(
            "multi_example expects a 48 x 48 actuator-space DM surface",
        ))
        dm_center = n_actuators / 2
        prop_dm(wfo, dm_surface, dm_center, dm_center, 1.0e-3)
    end

    prop_lens(wfo, fl_lens)
    prop_propagate(wfo, fl_lens)
    return prop_end(wfo; noabs=true)
end

if abspath(PROGRAM_FILE) == @__FILE__
    using Plots

    model = prepare_model(
        multi_example,
        0.55,
        256;
        name=:multi_example,
        PASSVALUE=[(; use_dm=false), (; use_dm=true, dm=zeros(48, 48))],
        pool_size=2,
    )
    stack, samplings = prop_run_multi(model)
    println("multi_example: samplings = ", samplings)
    intensity = abs2.(@view stack[:, :, 1])
    log_intensity = log10.(max.(intensity, eps(eltype(intensity))))
    display(heatmap(log_intensity; aspect_ratio=:equal, color=:grays, title="multi_example [0]"))
end
