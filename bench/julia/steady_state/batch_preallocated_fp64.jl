using BenchmarkTools
using FFTW
using JSON3
using Proper

include(joinpath(@__DIR__, "..", "..", "common", "metadata.jl"))
include(joinpath(@__DIR__, "..", "..", "common", "prepared_execution_workloads.jl"))
using .BenchMetadata

@inline function trial_stats(trial::BenchmarkTools.Trial)
    estimate = median(trial)
    return Dict(
        "median_ns" => estimate.time,
        "memory_bytes" => estimate.memory,
        "allocations" => estimate.allocs,
        "samples" => length(trial.times),
    )
end

function write_report(path::AbstractString, report)
    mkpath(dirname(path))
    open(path, "w") do io
        JSON3.write(io, report)
    end
    println(report)
    return report
end

grid_n = parse(Int, get(ENV, "PROPER_BENCH_GRID_N", string(PREPARED_STEADY_GRID_N)))
samples = parse(Int, get(ENV, "PROPER_BENCH_SAMPLES", string(PREPARED_STEADY_SAMPLES)))
report_path = get(
    ENV,
    "PROPER_BENCH_REPORT",
    joinpath(@__DIR__, "..", "..", "reports", "julia_batch_preallocated_cpu_fp64.json"),
)

FFTW.set_num_threads(parse(Int, get(ENV, "PROPER_BENCH_FFTW_THREADS", "1")))

allocating_runs = prepare_steady_state_sweep(
    Float64,
    cpu_prepared_context,
    PREPARED_BATCH_WAVELENGTHS,
    grid_n;
    name_prefix=:batch_allocating_cpu,
)
stack, samplings, prepared_runs = prepare_preallocated_cpu_sweep(
    PREPARED_BATCH_WAVELENGTHS,
    grid_n,
)

allocating_reference, reference_samplings = prop_run_multi(allocating_runs)
prop_run_multi!(stack, samplings, prepared_runs)
isapprox(stack, allocating_reference; atol=5e-13, rtol=5e-13) ||
    error("preallocated batch output differs from allocating reference")
samplings == reference_samplings ||
    error("preallocated batch samplings differ from allocating reference")

allocating_workload = () -> (prop_run_multi(allocating_runs); nothing)
preallocated_workload = () -> (
    prop_run_multi!(stack, samplings, prepared_runs);
    nothing
)
allocating_trial = benchmark_prepared_trial(allocating_workload; samples=samples)
preallocated_trial = benchmark_prepared_trial(preallocated_workload; samples=samples)

allocating_stats = trial_stats(allocating_trial)
preallocated_stats = trial_stats(preallocated_trial)
report = Dict(
    "meta" => merge(
        benchmark_metadata(run_tag="batch_preallocated_cpu_fp64", backend=:cpu),
        Dict(
            "grid_n" => grid_n,
            "precision" => "Float64",
            "batch_size" => length(PREPARED_BATCH_WAVELENGTHS),
            "julia_threads" => Threads.nthreads(),
            "fftw_threads" => FFTW.get_num_threads(),
        ),
    ),
    "policy" => "warmed allocating prop_run_multi versus caller-owned prop_run_multi!; preparation and TTFx excluded; identical output required before timing",
    "allocating" => allocating_stats,
    "preallocated" => preallocated_stats,
    "speedup" => allocating_stats["median_ns"] / preallocated_stats["median_ns"],
    "allocation_reduction" => 1 - preallocated_stats["memory_bytes"] / allocating_stats["memory_bytes"],
)

write_report(report_path, report)
