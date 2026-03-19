using JSON3
using Proper
include(joinpath(@__DIR__, "..", "..", "common", "metadata.jl"))
include(joinpath(@__DIR__, "..", "..", "common", "prepared_execution_workloads.jl"))
using .BenchMetadata

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

function write_cpu_report(path::AbstractString, report)
    mkpath(dirname(path))
    open(path, "w") do io
        JSON3.write(io, report)
    end
    println(report)
    return report
end

const RUN_TAG = "steady_state_fp32"
const REPORT_PATH = joinpath(@__DIR__, "..", "..", "reports", "julia_steady_state_fp32.json")

prepared = prepare_steady_state_model(Float32, PREPARED_STEADY_GRID_N, cpu_prepared_context(Float32, PREPARED_STEADY_GRID_N); name=:steady_state_cpu_fp32)
trial = benchmark_prepared_trial(() -> run_prepared_workload(prepared))
stats = trial_stats(trial)

report = Dict(
    "meta" => merge(benchmark_metadata(run_tag=RUN_TAG, backend=:cpu), Dict("grid_n" => PREPARED_STEADY_GRID_N, "precision" => "Float32")),
    "policy" => "steady-state CPU Float32 workload timing via BenchmarkTools with evals=1; TTFx excluded",
    "stats" => stats,
)

write_cpu_report(REPORT_PATH, report)
