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

const RUN_TAG = "batch_throughput_cpu"
const RUN_TAG_FP64 = "batch_throughput_cpu_fp64"
const REPORT_PATH = joinpath(@__DIR__, "..", "..", "reports", "julia_batch_throughput_cpu.json")
const REPORT_PATH_FP64 = joinpath(@__DIR__, "..", "..", "reports", "julia_batch_throughput_cpu_fp64.json")

prepared = prepare_steady_state_sweep(Float64, cpu_prepared_context, PREPARED_BATCH_WAVELENGTHS, PREPARED_STEADY_GRID_N; name_prefix=:batch_cpu)
trial = benchmark_prepared_trial(() -> run_prepared_sweep_workload(prepared))
stats = trial_stats(trial)

report = Dict(
    "meta" => merge(benchmark_metadata(run_tag=RUN_TAG, backend=:cpu), Dict("grid_n" => PREPARED_STEADY_GRID_N, "precision" => "Float64", "batch_size" => length(PREPARED_BATCH_WAVELENGTHS))),
    "policy" => "prepared batch throughput timing for a fixed wavelength sweep via prop_run_multi(vector_of_prepared_runs); TTFx excluded",
    "stats" => stats,
)

write_cpu_report(REPORT_PATH, report)
report_fp64 = copy(report)
meta_fp64 = copy(report["meta"])
meta_fp64["run_tag"] = RUN_TAG_FP64
report_fp64["meta"] = meta_fp64
write_cpu_report(REPORT_PATH_FP64, report_fp64)
