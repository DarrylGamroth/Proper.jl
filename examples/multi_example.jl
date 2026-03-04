using proper
using Plots

function prescription(lambda_m, n, passvalue; kwargs...)
    wf = prop_begin(2.4, lambda_m, n; beam_diam_fraction=0.5)
    prop_circular_aperture(wf, 0.6)
    prop_lens(wf, 10.0 + passvalue)
    prop_propagate(wf, 10.0)
    return wf
end

stack, samplings = prop_run_multi(prescription, 0.55, 256; PASSVALUE=[0.0, 1.0, 2.0])
println("multi_example: samplings = ", samplings)
heatmap(log10.(stack[:, :, 1] .+ eps()), aspect_ratio=:equal, title="multi_example [0]")
