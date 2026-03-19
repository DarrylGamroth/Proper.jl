using JSON3
include(joinpath(@__DIR__, "..", "..", "common", "metadata.jl"))
include(joinpath(@__DIR__, "..", "..", "common", "amdgpu_support.jl"))
using .BenchMetadata

const SUPPORTED_RUN_TAG = "steady_state_amdgpu_supported_kernels"
const SUPPORTED_REPORT_PATH = joinpath(@__DIR__, "..", "..", "reports", "amdgpu_supported_kernels.json")

reason = get(ENV, "AMDGPU_SKIP_REASON", "AMDGPU supported-kernel benchmark skipped")
write_benchmark_report(SUPPORTED_REPORT_PATH, skipped_amdgpu_report(SUPPORTED_RUN_TAG, reason))
