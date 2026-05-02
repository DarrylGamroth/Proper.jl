using BenchmarkTools
using Statistics

struct AMDGPUBenchmarkCase{SetupF,RunF}
    setup!::SetupF
    run_async!::RunF
end

const AMDGPU_WAVEFRONT_KERNEL_ORDER = (
    "prop_qphase",
    "prop_ptp",
    "prop_wts",
    "prop_stw",
    "prop_circular_aperture",
    "prop_dm_direct_map",
    "prop_end_mutating",
)

function run_amdgpu_benchmark_case(case::AMDGPUBenchmarkCase, samples::Integer)
    setup! = case.setup!
    run_async! = case.run_async!
    return run(@benchmarkable begin
        $run_async!()
        amdgpu_sync()
    end setup=($setup!()) evals=1 samples=samples)
end

function warmup_amdgpu_benchmark_case(case::AMDGPUBenchmarkCase, warmup_iters::Integer)
    for _ in 1:warmup_iters
        case.setup!()
        case.run_async!()
        amdgpu_sync()
    end
    return nothing
end

function build_amdgpu_wavefront_kernel_cases(::Type{T}, grid_n::Integer) where {T<:AbstractFloat}
    phase = T(10)
    distance = T(0.01)
    radius = T(0.6)
    diam_frac = T(0.5)

    wf_q = amdgpu_wavefront_begin(T, 2.4, 0.55e-6, grid_n; beam_diam_fraction=diam_frac)
    ctx_q = RunContext(wf_q)

    wf_p = amdgpu_wavefront_begin(T, 2.4, 0.55e-6, grid_n; beam_diam_fraction=diam_frac)
    ctx_p = RunContext(wf_p)

    wf_w = amdgpu_wavefront_begin(T, 2.4, 0.55e-6, grid_n; beam_diam_fraction=diam_frac)
    ctx_w = RunContext(wf_w)

    wf_s = amdgpu_wavefront_begin(T, 2.4, 0.55e-6, grid_n; beam_diam_fraction=diam_frac)
    ctx_s = RunContext(wf_s)
    prop_wts(wf_s, distance, ctx_s)

    wf_a = amdgpu_wavefront_begin(T, 2.4, 0.55e-6, grid_n; beam_diam_fraction=diam_frac)

    wf_dm = amdgpu_wavefront_begin(T, 2.4, 0.55e-6, grid_n; beam_diam_fraction=diam_frac)
    dm_map = AMDGPU.fill(T(1e-9), grid_n, grid_n)

    wf_e = amdgpu_wavefront_begin(T, 2.4, 0.55e-6, grid_n; beam_diam_fraction=diam_frac)
    prop_circular_aperture(wf_e, radius)
    out_end = amdgpu_zeros(T, grid_n, grid_n)

    snap_q = capture_wavefront_state(wf_q)
    snap_p = capture_wavefront_state(wf_p)
    snap_w = capture_wavefront_state(wf_w)
    snap_s = capture_wavefront_state(wf_s)
    snap_a = capture_wavefront_state(wf_a)
    snap_dm = capture_wavefront_state(wf_dm)
    snap_e = capture_wavefront_state(wf_e)
    amdgpu_sync()

    prop_qphase(wf_q, phase, ctx_q)
    wf_p.reference_surface = Proper.PLANAR
    prop_ptp(wf_p, distance, ctx_p)
    prop_wts(wf_w, distance, ctx_w)
    prop_stw(wf_s, distance, ctx_s)
    prop_circular_aperture(wf_a, radius)
    prop_dm(wf_dm, dm_map)
    prop_end!(out_end, wf_e)
    amdgpu_sync()

    return (
        prop_qphase=AMDGPUBenchmarkCase(
            () -> (restore_wavefront_state!(wf_q, snap_q); amdgpu_sync()),
            () -> prop_qphase(wf_q, phase, ctx_q),
        ),
        prop_ptp=AMDGPUBenchmarkCase(
            () -> (restore_wavefront_state!(wf_p, snap_p); amdgpu_sync()),
            () -> prop_ptp(wf_p, distance, ctx_p),
        ),
        prop_wts=AMDGPUBenchmarkCase(
            () -> (restore_wavefront_state!(wf_w, snap_w); amdgpu_sync()),
            () -> prop_wts(wf_w, distance, ctx_w),
        ),
        prop_stw=AMDGPUBenchmarkCase(
            () -> (restore_wavefront_state!(wf_s, snap_s); amdgpu_sync()),
            () -> prop_stw(wf_s, distance, ctx_s),
        ),
        prop_circular_aperture=AMDGPUBenchmarkCase(
            () -> (restore_wavefront_state!(wf_a, snap_a); amdgpu_sync()),
            () -> prop_circular_aperture(wf_a, radius),
        ),
        prop_dm_direct_map=AMDGPUBenchmarkCase(
            () -> (restore_wavefront_state!(wf_dm, snap_dm); amdgpu_sync()),
            () -> prop_dm(wf_dm, dm_map),
        ),
        prop_end_mutating=AMDGPUBenchmarkCase(
            () -> (restore_wavefront_state!(wf_e, snap_e); amdgpu_sync()),
            () -> prop_end!(out_end, wf_e),
        ),
    )
end

function benchmark_amdgpu_wavefront_kernel_stats(::Type{T}; grid_n::Integer, samples::Integer) where {T<:AbstractFloat}
    cases = build_amdgpu_wavefront_kernel_cases(T, grid_n)
    stats = Dict{String,Any}()
    for name in AMDGPU_WAVEFRONT_KERNEL_ORDER
        case = getproperty(cases, Symbol(name))
        warmup_amdgpu_benchmark_case(case, 2)
        stats[name] = trial_stats(run_amdgpu_benchmark_case(case, samples))
    end
    return stats
end
