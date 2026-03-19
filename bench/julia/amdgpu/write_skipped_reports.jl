using JSON3
include(joinpath(@__DIR__, "..", "..", "common", "metadata.jl"))
include(joinpath(@__DIR__, "..", "..", "common", "amdgpu_support.jl"))
using .BenchMetadata

const STEADY_RUN_TAG = "steady_state_amdgpu"
const STEADY_FP64_RUN_TAG = "steady_state_amdgpu_fp64"
const STEADY_FP32_RUN_TAG = "steady_state_amdgpu_fp32"
const SUPPORTED_RUN_TAG = "steady_state_amdgpu_supported_kernels"
const CORE_TAIL_RUN_TAG = "core_propagation_tail_amdgpu"
const BATCH_RUN_TAG = "batch_throughput_amdgpu"
const BATCH_FP64_RUN_TAG = "batch_throughput_amdgpu_fp64"
const BATCH_FP32_RUN_TAG = "batch_throughput_amdgpu_fp32"
const STEADY_REPORT_PATH = joinpath(@__DIR__, "..", "..", "reports", "julia_amdgpu_steady_state.json")
const STEADY_FP64_REPORT_PATH = joinpath(@__DIR__, "..", "..", "reports", "julia_amdgpu_steady_state_fp64.json")
const STEADY_FP32_REPORT_PATH = joinpath(@__DIR__, "..", "..", "reports", "julia_amdgpu_steady_state_fp32.json")
const SUPPORTED_REPORT_PATH = joinpath(@__DIR__, "..", "..", "reports", "amdgpu_supported_kernels.json")
const CORE_TAIL_REPORT_PATH = joinpath(@__DIR__, "..", "..", "reports", "julia_core_propagation_tail_amdgpu.json")
const BATCH_REPORT_PATH = joinpath(@__DIR__, "..", "..", "reports", "julia_batch_throughput_amdgpu.json")
const BATCH_FP64_REPORT_PATH = joinpath(@__DIR__, "..", "..", "reports", "julia_batch_throughput_amdgpu_fp64.json")
const BATCH_FP32_REPORT_PATH = joinpath(@__DIR__, "..", "..", "reports", "julia_batch_throughput_amdgpu_fp32.json")

reason = get(ENV, "AMDGPU_SKIP_REASON", "AMDGPU benchmark skipped")
write_benchmark_report(STEADY_REPORT_PATH, skipped_amdgpu_report(STEADY_RUN_TAG, reason))
write_benchmark_report(STEADY_FP64_REPORT_PATH, skipped_amdgpu_report(STEADY_FP64_RUN_TAG, reason))
write_benchmark_report(STEADY_FP32_REPORT_PATH, skipped_amdgpu_report(STEADY_FP32_RUN_TAG, reason))
write_benchmark_report(SUPPORTED_REPORT_PATH, skipped_amdgpu_report(SUPPORTED_RUN_TAG, reason))
write_benchmark_report(CORE_TAIL_REPORT_PATH, skipped_amdgpu_report(CORE_TAIL_RUN_TAG, reason))
write_benchmark_report(BATCH_REPORT_PATH, skipped_amdgpu_report(BATCH_RUN_TAG, reason))
write_benchmark_report(BATCH_FP64_REPORT_PATH, skipped_amdgpu_report(BATCH_FP64_RUN_TAG, reason))
write_benchmark_report(BATCH_FP32_REPORT_PATH, skipped_amdgpu_report(BATCH_FP32_RUN_TAG, reason))
