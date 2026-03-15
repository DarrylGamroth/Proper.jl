using BenchmarkTools
using JSON3
using Proper
include(joinpath(@__DIR__, "..", "..", "common", "workloads.jl"))
include(joinpath(@__DIR__, "..", "..", "common", "metadata.jl"))
using .Workloads
using .BenchMetadata

function steady_state_prescription(λm, n; kwargs...)
    wf = prop_begin(2.4, λm, n; beam_diam_fraction=0.5)
    prop_circular_aperture(wf, 0.6)
    prop_lens(wf, 20.0)
    prop_propagate(wf, 20.0)
    return wf
end

const PREPARED_STEADY_STATE = prepare_prescription(steady_state_prescription, 0.55, 512)

workload() = prop_run(PREPARED_STEADY_STATE)

# Warmup: exclude compilation from steady-state timings.
workload()

b = @benchmark workload()
stats = Dict(
    "median_ns" => median(b).time,
    "mean_ns" => mean(b).time,
    "min_ns" => minimum(b).time,
    "max_ns" => maximum(b).time,
    "samples" => length(b.times),
)

report = Dict(
    "meta" => benchmark_metadata(run_tag="steady_state"),
    "policy" => "TTFx excluded from these timings",
    "stats" => stats,
)

mkpath(joinpath(@__DIR__, "..", "..", "reports"))
out = joinpath(@__DIR__, "..", "..", "reports", "julia_steady_state.json")
open(out, "w") do io
    JSON3.write(io, report)
end
println(report)
