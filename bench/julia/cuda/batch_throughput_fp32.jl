using JSON3
using Proper
using CUDA
include(joinpath(@__DIR__, "..", "..", "common", "metadata.jl"))
include(joinpath(@__DIR__, "..", "..", "common", "cuda_support.jl"))
include(joinpath(@__DIR__, "..", "..", "common", "prepared_execution_workloads.jl"))
using .BenchMetadata

const RUN_TAG = "batch_throughput_cuda_fp32"
const REPORT_PATH = joinpath(@__DIR__, "..", "..", "reports", "julia_batch_throughput_cuda_fp32.json")

CUDA.allowscalar(false)
prepared = prepare_steady_state_sweep(Float32, (T, n) -> backend_prepared_context(cuda_wavefront_begin, T, n), PREPARED_BATCH_WAVELENGTHS, PREPARED_STEADY_GRID_N; name_prefix=:batch_cuda_fp32)
trial = benchmark_prepared_trial(() -> run_prepared_sweep_workload(prepared, cuda_sync))

report = Dict(
    "meta" => merge(cuda_report_meta(RUN_TAG; device=cuda_device_label()), Dict("grid_n" => PREPARED_STEADY_GRID_N, "precision" => "Float32", "batch_size" => length(PREPARED_BATCH_WAVELENGTHS))),
    "policy" => "prepared CUDA Float32 batch throughput timing for a fixed wavelength sweep via prop_run_multi(vector_of_prepared_runs); TTFx and initial CUDA context setup excluded; per-sample synchronization included",
    "stats" => trial_stats(trial),
)

write_benchmark_report(REPORT_PATH, report)
