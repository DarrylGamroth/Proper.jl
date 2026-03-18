using BenchmarkTools
using JSON3
using Proper
include(joinpath(@__DIR__, "..", "..", "common", "metadata.jl"))
include(joinpath(@__DIR__, "..", "..", "common", "wavefront_state.jl"))
include(joinpath(@__DIR__, "..", "..", "common", "core_propagation_tail.jl"))
using .BenchMetadata

const RUN_TAG = "core_propagation_tail_cpu"
const REPORT_PATH = joinpath(@__DIR__, "..", "..", "reports", "julia_core_propagation_tail_cpu.json")

@inline function trial_stats(t::BenchmarkTools.Trial)
    est = median(t)
    return Dict(
        "median_ns" => est.time,
        "median_allocs" => est.allocs,
        "median_bytes" => est.memory,
        "samples" => length(t.times),
    )
end

wf, ctx, snap = cpu_core_propagation_tail_case(Float64)
restore_and_run_core_propagation_tail!(wf, snap, ctx)

trial = run(@benchmarkable begin
    restore_and_run_core_propagation_tail!($wf, $snap, $ctx)
end evals=1 samples=CORE_PROPAGATION_TAIL_SAMPLES)

report = Dict(
    "meta" => merge(
        benchmark_metadata(run_tag=RUN_TAG, backend=:cpu),
        Dict("grid_n" => CORE_PROPAGATION_TAIL_GRID_N, "precision" => "Float64"),
    ),
    "policy" => "steady-state synthetic core propagation tail timing via BenchmarkTools with evals=1 and per-sample state restore; TTFx excluded",
    "sequence" => [Dict("name" => name, "op" => String(op), "value" => value) for (name, op, value) in CORE_PROPAGATION_TAIL_SEQUENCE],
    "stats" => trial_stats(trial),
)

mkpath(joinpath(@__DIR__, "..", "..", "reports"))
open(REPORT_PATH, "w") do io
    JSON3.write(io, report)
end
println(report)
