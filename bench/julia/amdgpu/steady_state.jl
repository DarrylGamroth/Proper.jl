using JSON3
using Proper
using AMDGPU
include(joinpath(@__DIR__, "..", "..", "common", "metadata.jl"))
include(joinpath(@__DIR__, "..", "..", "common", "amdgpu_support.jl"))
include(joinpath(@__DIR__, "..", "..", "common", "amdgpu_steady_state_workload.jl"))
using .BenchMetadata

const RUN_TAG = "steady_state_amdgpu"
const RUN_TAG_FP64 = "steady_state_amdgpu_fp64"
const REPORT_PATH = joinpath(@__DIR__, "..", "..", "reports", "julia_amdgpu_steady_state.json")
const REPORT_PATH_FP64 = joinpath(@__DIR__, "..", "..", "reports", "julia_amdgpu_steady_state_fp64.json")

AMDGPU.allowscalar(false)
run_amdgpu_steady_state_report(Float64, RUN_TAG, REPORT_PATH; alias_specs=(REPORT_PATH_FP64 => RUN_TAG_FP64,))
