using Profile
using Proper

const REPORT_ROOT = joinpath(@__DIR__, "..", "..", "reports")
const BACKEND = isempty(ARGS) ? "cpu" : first(ARGS)

include(joinpath(@__DIR__, "..", "..", "common", "wavefront_state.jl"))
include(joinpath(@__DIR__, "..", "..", "common", "core_propagation_tail.jl"))

if BACKEND == "cuda"
    @eval using CUDA
    include(joinpath(@__DIR__, "..", "..", "common", "cuda_support.jl"))
    CUDA.allowscalar(false)
elseif BACKEND == "amdgpu"
    @eval using AMDGPU
    include(joinpath(@__DIR__, "..", "..", "common", "amdgpu_support.jl"))
    AMDGPU.allowscalar(false)
end

function backend_case(backend::String)
    if backend == "cpu"
        wf, ctx, snap = cpu_core_propagation_tail_case(Float64)
        sync! = () -> nothing
        label = "cpu"
    elseif backend == "cuda"
        wf, ctx, snap = gpu_core_propagation_tail_case(cuda_wavefront_begin, cuda_sync, Float64)
        sync! = cuda_sync
        label = "cuda"
    elseif backend == "amdgpu"
        wf, ctx, snap = gpu_core_propagation_tail_case(amdgpu_wavefront_begin, amdgpu_sync, Float64)
        sync! = amdgpu_sync
        label = "amdgpu"
    else
        throw(ArgumentError("backend must be one of: cpu, cuda, amdgpu"))
    end
    return wf, ctx, snap, sync!, label
end

function main(args)
    backend = isempty(args) ? "cpu" : first(args)
    iters = length(args) >= 2 ? parse(Int, args[2]) : CORE_PROPAGATION_TAIL_PROFILE_ITERS
    wf, ctx, snap, sync!, label = backend_case(backend)

    restore_and_run_core_propagation_tail!(wf, snap, ctx)
    sync!()

    Profile.clear()
    Profile.init()
    Profile.@profile begin
        for _ in 1:iters
            restore_and_run_core_propagation_tail!(wf, snap, ctx)
            sync!()
        end
    end

    mkpath(REPORT_ROOT)
    out = joinpath(REPORT_ROOT, "core_propagation_tail_profile_" * label * ".txt")
    open(out, "w") do io
        println(io, "Core propagation tail profile")
        println(io, "backend: ", label)
        println(io, "iterations: ", iters)
        println(io, "grid_n: ", CORE_PROPAGATION_TAIL_GRID_N)
        println(io, "sequence:")
        for (name, op, value) in CORE_PROPAGATION_TAIL_SEQUENCE
            println(io, "  - ", name, ": ", op, "(", value, ")")
        end
        println(io)
        Profile.print(io; sortedby=:count, maxdepth=24)
    end
    println(out)
end

main(ARGS)
