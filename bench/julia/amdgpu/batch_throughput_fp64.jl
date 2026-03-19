using JSON3
using Proper
using AMDGPU
include(joinpath(@__DIR__, "..", "..", "common", "metadata.jl"))
include(joinpath(@__DIR__, "..", "..", "common", "amdgpu_support.jl"))
include(joinpath(@__DIR__, "..", "..", "common", "prepared_execution_workloads.jl"))
using .BenchMetadata

const RUN_TAG = "batch_throughput_amdgpu"
const RUN_TAG_FP64 = "batch_throughput_amdgpu_fp64"
const REPORT_PATH = joinpath(@__DIR__, "..", "..", "reports", "julia_batch_throughput_amdgpu.json")
const REPORT_PATH_FP64 = joinpath(@__DIR__, "..", "..", "reports", "julia_batch_throughput_amdgpu_fp64.json")

AMDGPU.allowscalar(false)
prepared = prepare_steady_state_sweep(Float64, (T, n) -> backend_prepared_context(amdgpu_wavefront_begin, T, n), PREPARED_BATCH_WAVELENGTHS, PREPARED_STEADY_GRID_N; name_prefix=:batch_amdgpu)
trial = benchmark_prepared_trial(() -> run_prepared_sweep_workload(prepared, amdgpu_sync))

report = Dict(
    "meta" => merge(amdgpu_report_meta(RUN_TAG; device=amdgpu_device_label()), Dict("grid_n" => PREPARED_STEADY_GRID_N, "precision" => "Float64", "batch_size" => length(PREPARED_BATCH_WAVELENGTHS))),
    "policy" => "prepared AMDGPU batch throughput timing for a fixed wavelength sweep via prop_run_multi(vector_of_prepared_runs); TTFx and initial GPU context setup excluded; per-sample synchronization included",
    "stats" => trial_stats(trial),
)

write_benchmark_report(REPORT_PATH, report)
report_fp64 = copy(report)
meta_fp64 = copy(report["meta"])
meta_fp64["run_tag"] = RUN_TAG_FP64
report_fp64["meta"] = meta_fp64
write_benchmark_report(REPORT_PATH_FP64, report_fp64)
