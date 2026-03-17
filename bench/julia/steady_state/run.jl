using BenchmarkTools
using JSON3
using Proper
include(joinpath(@__DIR__, "..", "..", "common", "workloads.jl"))
include(joinpath(@__DIR__, "..", "..", "common", "metadata.jl"))
using .Workloads
using .BenchMetadata

const SAMPLES = 20

@inline function trial_stats(t::BenchmarkTools.Trial)
    est_med = median(t)
    est_mean = mean(t)
    est_min = minimum(t)
    est_max = maximum(t)
    return Dict(
        "median_ns" => est_med.time,
        "mean_ns" => est_mean.time,
        "min_ns" => est_min.time,
        "max_ns" => est_max.time,
        "samples" => length(t.times),
    )
end

function steady_state_prescription(λm, n; kwargs...)
    wf = prop_begin(2.4, λm, n; beam_diam_fraction=0.5)
    prop_circular_aperture(wf, 0.6)
    prop_lens(wf, 20.0)
    prop_propagate(wf, 20.0)
    return wf
end

const PREPARED_STEADY_STATE = prepare_model(:steady_state_cpu, steady_state_prescription, 0.55, 512; pool_size=1)

workload() = prop_run(PREPARED_STEADY_STATE)

# Warmup: exclude compilation from steady-state timings.
workload()

b = run(@benchmarkable workload() evals=1 samples=SAMPLES)
stats = trial_stats(b)

report = Dict(
    "meta" => benchmark_metadata(run_tag="steady_state"),
    "policy" => "steady-state CPU workload timing via BenchmarkTools with evals=1; TTFx excluded",
    "stats" => stats,
)

mkpath(joinpath(@__DIR__, "..", "..", "reports"))
out = joinpath(@__DIR__, "..", "..", "reports", "julia_steady_state.json")
open(out, "w") do io
    JSON3.write(io, report)
end
println(report)
