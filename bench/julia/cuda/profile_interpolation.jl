using Profile
using Proper
using CUDA
include(joinpath(@__DIR__, "..", "..", "common", "cuda_support.jl"))

const GRID_N = 512
const HOST_PROFILE_ITERS = 20
const DEVICE_PROFILE_ITERS = 100
const CASE_ORDER = ("cubic_conv_grid", "rotate", "resamplemap", "magnify_quick")

CUDA.allowscalar(false)

@inline function _case_names(args)
    isempty(args) && return CASE_ORDER
    requested = Tuple(args)
    invalid = filter(name -> !(name in CASE_ORDER), requested)
    isempty(invalid) || throw(ArgumentError("unknown case(s): $(join(invalid, ", "))"))
    return requested
end

function build_cubic_conv_case()
    img = CUDA.rand(Float32, GRID_N, GRID_N)
    out = similar(img)
    x = CUDA.CuArray(Float32.(1:GRID_N))
    y = CUDA.CuArray(Float32.(1:GRID_N))
    action() = Proper.prop_cubic_conv_grid!(out, img, x, y)
    return action
end

function build_rotate_case()
    img = CUDA.rand(Float32, GRID_N, GRID_N)
    out = similar(img)
    ctx = RunContext(typeof(img))
    action() = Proper.prop_rotate!(out, img, 5.0f0, ctx)
    return action
end

function build_resamplemap_case()
    dmap = CUDA.rand(Float32, GRID_N, GRID_N)
    wf = Proper.WaveFront(CUDA.fill(ComplexF32(1), GRID_N, GRID_N), 500f-9, 1f-3, 0f0, 1f0)
    ctx = RunContext(typeof(dmap))
    opts = Proper.ResampleMapOptions(wf, wf.sampling_m, Float32(GRID_N ÷ 2), Float32(GRID_N ÷ 2))
    out = similar(dmap, Float64, size(wf.field)...)
    action() = Proper.prop_resamplemap!(out, wf, dmap, opts, ctx)
    return action
end

function build_magnify_quick_case()
    img = CUDA.rand(Float32, GRID_N, GRID_N)
    out = similar(img)
    ctx = RunContext(typeof(img))
    action() = Proper.prop_magnify!(out, img, 1.15f0, ctx; QUICK=true)
    return action
end

function build_case(name::AbstractString)
    if name == "cubic_conv_grid"
        return build_cubic_conv_case()
    elseif name == "rotate"
        return build_rotate_case()
    elseif name == "resamplemap"
        return build_resamplemap_case()
    elseif name == "magnify_quick"
        return build_magnify_quick_case()
    end
    throw(ArgumentError("unknown case: $name"))
end

function warmup!(action)
    action()
    cuda_sync()
    return nothing
end

function host_profile!(name::AbstractString, action)
    println("== Host Profile: ", name, " ==")
    Profile.clear()
    Profile.@profile begin
        for _ in 1:HOST_PROFILE_ITERS
            action()
            cuda_sync()
        end
    end
    Profile.print()
    println()
    return nothing
end

function device_profile!(name::AbstractString, action)
    println("== Device Profile Marker: ", name, " ==")
    println("Use the CUDA profiler output below to decide whether tiling/shared-memory work is justified.")
    CUDA.@profile begin
        for _ in 1:DEVICE_PROFILE_ITERS
            action()
        end
        cuda_sync()
    end
    println()
    return nothing
end

function main(args=ARGS)
    println("CUDA interpolation profiling harness")
    println("grid_n = ", GRID_N)
    println("host_profile_iters = ", HOST_PROFILE_ITERS)
    println("device_profile_iters = ", DEVICE_PROFILE_ITERS)
    println()
    println("Interpretation:")
    println("- Use host profile output to confirm wrapper/dispatch/allocation overhead is not dominant.")
    println("- Use CUDA profiler output to decide whether interpolation kernels show enough memory reuse pressure to justify tiling/@localmem.")
    println("- Do not start tiling work for propagation FFT kernels from this harness; it is only for interpolation-family kernels.")
    println()

    for name in _case_names(args)
        action = build_case(name)
        warmup!(action)
        host_profile!(name, action)
        device_profile!(name, action)
    end

    return nothing
end

main()
