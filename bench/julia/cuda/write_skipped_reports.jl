using JSON3
include(joinpath(@__DIR__, "..", "..", "common", "metadata.jl"))
include(joinpath(@__DIR__, "..", "..", "common", "cuda_support.jl"))
using .BenchMetadata

const STEADY_RUN_TAG = "steady_state_cuda"
const STEADY_FP64_RUN_TAG = "steady_state_cuda_fp64"
const STEADY_FP32_RUN_TAG = "steady_state_cuda_fp32"
const SUPPORTED_RUN_TAG = "steady_state_cuda_supported_kernels"
const PRECISION_RUN_TAG = "steady_state_cuda_precision_split"
const ISOLATED_RUN_TAG = "cuda_isolated_wavefront_kernels"
const CORE_TAIL_RUN_TAG = "core_propagation_tail_cuda"
const BATCH_RUN_TAG = "batch_throughput_cuda"
const BATCH_FP64_RUN_TAG = "batch_throughput_cuda_fp64"
const BATCH_FP32_RUN_TAG = "batch_throughput_cuda_fp32"
const STEADY_REPORT_PATH = joinpath(@__DIR__, "..", "..", "reports", "julia_cuda_steady_state.json")
const STEADY_FP64_REPORT_PATH = joinpath(@__DIR__, "..", "..", "reports", "julia_cuda_steady_state_fp64.json")
const STEADY_FP32_REPORT_PATH = joinpath(@__DIR__, "..", "..", "reports", "julia_cuda_steady_state_fp32.json")
const SUPPORTED_REPORT_PATH = joinpath(@__DIR__, "..", "..", "reports", "cuda_supported_kernels.json")
const PRECISION_REPORT_PATH = joinpath(@__DIR__, "..", "..", "reports", "cuda_precision_split.json")
const ISOLATED_REPORT_PATH = joinpath(@__DIR__, "..", "..", "reports", "cuda_isolated_wavefront_kernels.json")
const CORE_TAIL_REPORT_PATH = joinpath(@__DIR__, "..", "..", "reports", "julia_core_propagation_tail_cuda.json")
const BATCH_REPORT_PATH = joinpath(@__DIR__, "..", "..", "reports", "julia_batch_throughput_cuda.json")
const BATCH_FP64_REPORT_PATH = joinpath(@__DIR__, "..", "..", "reports", "julia_batch_throughput_cuda_fp64.json")
const BATCH_FP32_REPORT_PATH = joinpath(@__DIR__, "..", "..", "reports", "julia_batch_throughput_cuda_fp32.json")

reason = get(ENV, "CUDA_SKIP_REASON", "CUDA benchmark skipped")
write_benchmark_report(STEADY_REPORT_PATH, skipped_cuda_report(STEADY_RUN_TAG, reason))
write_benchmark_report(STEADY_FP64_REPORT_PATH, skipped_cuda_report(STEADY_FP64_RUN_TAG, reason))
write_benchmark_report(STEADY_FP32_REPORT_PATH, skipped_cuda_report(STEADY_FP32_RUN_TAG, reason))
write_benchmark_report(SUPPORTED_REPORT_PATH, skipped_cuda_report(SUPPORTED_RUN_TAG, reason))
write_benchmark_report(PRECISION_REPORT_PATH, skipped_cuda_report(PRECISION_RUN_TAG, reason))
write_benchmark_report(ISOLATED_REPORT_PATH, skipped_cuda_report(ISOLATED_RUN_TAG, reason))
write_benchmark_report(CORE_TAIL_REPORT_PATH, skipped_cuda_report(CORE_TAIL_RUN_TAG, reason))
write_benchmark_report(BATCH_REPORT_PATH, skipped_cuda_report(BATCH_RUN_TAG, reason))
write_benchmark_report(BATCH_FP64_REPORT_PATH, skipped_cuda_report(BATCH_FP64_RUN_TAG, reason))
write_benchmark_report(BATCH_FP32_REPORT_PATH, skipped_cuda_report(BATCH_FP32_RUN_TAG, reason))
