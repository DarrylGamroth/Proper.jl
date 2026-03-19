using JSON3
using Proper
using AMDGPU
include(joinpath(@__DIR__, "..", "..", "common", "metadata.jl"))
include(joinpath(@__DIR__, "..", "..", "common", "amdgpu_support.jl"))
include(joinpath(@__DIR__, "..", "..", "common", "amdgpu_steady_state_workload.jl"))
using .BenchMetadata

const RUN_TAG = "steady_state_amdgpu_fp32"
const REPORT_PATH = joinpath(@__DIR__, "..", "..", "reports", "julia_amdgpu_steady_state_fp32.json")

AMDGPU.allowscalar(false)
run_amdgpu_steady_state_report(Float32, RUN_TAG, REPORT_PATH)
