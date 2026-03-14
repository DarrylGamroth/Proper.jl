using BenchmarkTools
using CUDA
using Statistics

struct CUDABenchmarkCase{SetupF,RunF}
    setup!::SetupF
    run_async!::RunF
end

const CUDA_WAVEFRONT_KERNEL_ORDER = (
    "prop_qphase",
    "prop_ptp",
    "prop_wts",
    "prop_stw",
    "prop_circular_aperture",
    "prop_end_mutating",
)

function run_cuda_benchmark_case(case::CUDABenchmarkCase, samples::Integer)
    setup! = case.setup!
    run_async! = case.run_async!
    return run(@benchmarkable begin
        $run_async!()
        cuda_sync()
    end setup=($setup!()) evals=1 samples=samples)
end

function warmup_cuda_benchmark_case(case::CUDABenchmarkCase, warmup_iters::Integer)
    for _ in 1:warmup_iters
        case.setup!()
        case.run_async!()
        cuda_sync()
    end
    return nothing
end

function sample_stats(values)
    isempty(values) && return Dict("samples" => 0)
    data = Float64.(values)
    return Dict(
        "median_ns" => median(data),
        "min_ns" => minimum(data),
        "max_ns" => maximum(data),
        "samples" => length(data),
    )
end

function collect_cuda_host_device_samples(case::CUDABenchmarkCase, samples::Integer)
    host_ns = Vector{Float64}(undef, samples)
    device_ns = Vector{Float64}(undef, samples)

    for i in 1:samples
        case.setup!()
        cuda_sync()
        host_t0 = time_ns()
        device_elapsed_s = CUDA.@elapsed begin
            case.run_async!()
        end
        cuda_sync()
        host_ns[i] = time_ns() - host_t0
        device_ns[i] = device_elapsed_s * 1e9
    end

    return Dict(
        "host_wall" => sample_stats(host_ns),
        "device" => sample_stats(device_ns),
    )
end

function build_cuda_wavefront_kernel_cases(::Type{T}, grid_n::Integer) where {T<:AbstractFloat}
    phase = T(10)
    distance = T(0.01)
    radius = T(0.6)
    diam_frac = T(0.5)

    wf_q = cuda_wavefront_begin(T, 2.4, 0.55e-6, grid_n; beam_diam_fraction=diam_frac)
    ctx_q = RunContext(wf_q)

    wf_p = cuda_wavefront_begin(T, 2.4, 0.55e-6, grid_n; beam_diam_fraction=diam_frac)
    ctx_p = RunContext(wf_p)

    wf_w = cuda_wavefront_begin(T, 2.4, 0.55e-6, grid_n; beam_diam_fraction=diam_frac)
    ctx_w = RunContext(wf_w)

    wf_s = cuda_wavefront_begin(T, 2.4, 0.55e-6, grid_n; beam_diam_fraction=diam_frac)
    ctx_s = RunContext(wf_s)
    prop_wts(wf_s, distance, ctx_s)

    wf_a = cuda_wavefront_begin(T, 2.4, 0.55e-6, grid_n; beam_diam_fraction=diam_frac)

    wf_e = cuda_wavefront_begin(T, 2.4, 0.55e-6, grid_n; beam_diam_fraction=diam_frac)
    prop_circular_aperture(wf_e, radius)
    out_end = cuda_zeros(T, grid_n, grid_n)

    snap_q = capture_wavefront_state(wf_q)
    snap_p = capture_wavefront_state(wf_p)
    snap_w = capture_wavefront_state(wf_w)
    snap_s = capture_wavefront_state(wf_s)
    snap_a = capture_wavefront_state(wf_a)
    snap_e = capture_wavefront_state(wf_e)
    cuda_sync()

    # Warmup: exclude compilation and initial CUDA context setup.
    prop_qphase(wf_q, phase, ctx_q)
    wf_p.reference_surface = Proper.PLANAR
    prop_ptp(wf_p, distance, ctx_p)
    prop_wts(wf_w, distance, ctx_w)
    prop_stw(wf_s, distance, ctx_s)
    prop_circular_aperture(wf_a, radius)
    prop_end!(out_end, wf_e)
    cuda_sync()

    return (
        prop_qphase=CUDABenchmarkCase(
            () -> (restore_wavefront_state!(wf_q, snap_q); cuda_sync()),
            () -> prop_qphase(wf_q, phase, ctx_q),
        ),
        prop_ptp=CUDABenchmarkCase(
            () -> (restore_wavefront_state!(wf_p, snap_p); cuda_sync()),
            () -> prop_ptp(wf_p, distance, ctx_p),
        ),
        prop_wts=CUDABenchmarkCase(
            () -> (restore_wavefront_state!(wf_w, snap_w); cuda_sync()),
            () -> prop_wts(wf_w, distance, ctx_w),
        ),
        prop_stw=CUDABenchmarkCase(
            () -> (restore_wavefront_state!(wf_s, snap_s); cuda_sync()),
            () -> prop_stw(wf_s, distance, ctx_s),
        ),
        prop_circular_aperture=CUDABenchmarkCase(
            () -> (restore_wavefront_state!(wf_a, snap_a); cuda_sync()),
            () -> prop_circular_aperture(wf_a, radius),
        ),
        prop_end_mutating=CUDABenchmarkCase(
            () -> (restore_wavefront_state!(wf_e, snap_e); cuda_sync()),
            () -> prop_end!(out_end, wf_e),
        ),
    )
end

function benchmark_cuda_wavefront_kernel_stats(::Type{T}; grid_n::Integer, samples::Integer) where {T<:AbstractFloat}
    cases = build_cuda_wavefront_kernel_cases(T, grid_n)
    stats = Dict{String,Any}()
    for name in CUDA_WAVEFRONT_KERNEL_ORDER
        case = getproperty(cases, Symbol(name))
        stats[name] = trial_stats(run_cuda_benchmark_case(case, samples))
    end
    return stats
end
