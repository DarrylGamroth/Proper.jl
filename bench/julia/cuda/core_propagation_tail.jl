using BenchmarkTools
using JSON3
using Proper
using CUDA
include(joinpath(@__DIR__, "..", "..", "common", "metadata.jl"))
include(joinpath(@__DIR__, "..", "..", "common", "cuda_support.jl"))
include(joinpath(@__DIR__, "..", "..", "common", "wavefront_state.jl"))
include(joinpath(@__DIR__, "..", "..", "common", "core_propagation_tail.jl"))
using .BenchMetadata

const RUN_TAG = "core_propagation_tail_cuda"
const REPORT_PATH = joinpath(@__DIR__, "..", "..", "reports", "julia_core_propagation_tail_cuda.json")

CUDA.allowscalar(false)

wf, ctx, snap = gpu_core_propagation_tail_case(cuda_wavefront_begin, cuda_sync, Float64)
restore_and_run_core_propagation_tail!(wf, snap, ctx)
restore_and_run_core_propagation_tail!(wf, snap, ctx)
cuda_sync()

trial = run(@benchmarkable begin
    restore_and_run_core_propagation_tail!($wf, $snap, $ctx)
    cuda_sync()
end evals=1 samples=CORE_PROPAGATION_TAIL_SAMPLES)

report = Dict(
    "meta" => merge(
        cuda_report_meta(RUN_TAG; device=cuda_device_label()),
        Dict("grid_n" => CORE_PROPAGATION_TAIL_GRID_N, "precision" => "Float64"),
    ),
    "policy" => "steady-state synthetic core propagation tail timing via BenchmarkTools with evals=1, per-sample state restore, and per-sample synchronization; TTFx excluded",
    "sequence" => [Dict("name" => name, "op" => String(op), "value" => value) for (name, op, value) in CORE_PROPAGATION_TAIL_SEQUENCE],
    "stats" => trial_stats(trial),
)

write_benchmark_report(REPORT_PATH, report)
