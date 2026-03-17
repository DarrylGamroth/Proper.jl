using JSON3
include(joinpath(@__DIR__, "..", "..", "common", "metadata.jl"))
include(joinpath(@__DIR__, "..", "..", "common", "amdgpu_support.jl"))
using .BenchMetadata

const STEADY_RUN_TAG = "steady_state_amdgpu"
const SUPPORTED_RUN_TAG = "steady_state_amdgpu_supported_kernels"
const STEADY_REPORT_PATH = joinpath(@__DIR__, "..", "..", "reports", "julia_amdgpu_steady_state.json")
const SUPPORTED_REPORT_PATH = joinpath(@__DIR__, "..", "..", "reports", "amdgpu_supported_kernels.json")

reason = get(ENV, "AMDGPU_SKIP_REASON", "AMDGPU benchmark skipped")
write_benchmark_report(STEADY_REPORT_PATH, skipped_amdgpu_report(STEADY_RUN_TAG, reason))
write_benchmark_report(SUPPORTED_REPORT_PATH, skipped_amdgpu_report(SUPPORTED_RUN_TAG, reason))
