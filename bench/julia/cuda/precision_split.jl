using BenchmarkTools
using JSON3
using Proper
using CUDA
include(joinpath(@__DIR__, "..", "..", "common", "metadata.jl"))
include(joinpath(@__DIR__, "..", "..", "common", "cuda_support.jl"))
include(joinpath(@__DIR__, "..", "..", "common", "wavefront_state.jl"))
using .BenchMetadata

const RUN_TAG = "steady_state_cuda_precision_split"
const REPORT_PATH = joinpath(@__DIR__, "..", "..", "reports", "cuda_precision_split.json")
const GRID_N = 512
const N_SAMPLES = 20

CUDA.allowscalar(false)

function benchmark_precision_split()
    wf_q64 = cuda_wavefront_begin(Float64, 2.4, 0.55e-6, GRID_N; beam_diam_fraction=0.5)
    ctx_q64 = RunContext(wf_q64)
    wf_q32 = cuda_wavefront_begin(Float32, 2.4, 0.55e-6, GRID_N; beam_diam_fraction=0.5f0)
    ctx_q32 = RunContext(wf_q32)

    wf_p64 = cuda_wavefront_begin(Float64, 2.4, 0.55e-6, GRID_N; beam_diam_fraction=0.5)
    ctx_p64 = RunContext(wf_p64)
    wf_p32 = cuda_wavefront_begin(Float32, 2.4, 0.55e-6, GRID_N; beam_diam_fraction=0.5f0)
    ctx_p32 = RunContext(wf_p32)

    wf_w64 = cuda_wavefront_begin(Float64, 2.4, 0.55e-6, GRID_N; beam_diam_fraction=0.5)
    ctx_w64 = RunContext(wf_w64)
    wf_w32 = cuda_wavefront_begin(Float32, 2.4, 0.55e-6, GRID_N; beam_diam_fraction=0.5f0)
    ctx_w32 = RunContext(wf_w32)

    wf_s64 = cuda_wavefront_begin(Float64, 2.4, 0.55e-6, GRID_N; beam_diam_fraction=0.5)
    ctx_s64 = RunContext(wf_s64)
    wf_s32 = cuda_wavefront_begin(Float32, 2.4, 0.55e-6, GRID_N; beam_diam_fraction=0.5f0)
    ctx_s32 = RunContext(wf_s32)
    prop_wts(wf_s64, 0.01, ctx_s64)
    prop_wts(wf_s32, 0.01f0, ctx_s32)

    wf_a64 = cuda_wavefront_begin(Float64, 2.4, 0.55e-6, GRID_N; beam_diam_fraction=0.5)
    wf_a32 = cuda_wavefront_begin(Float32, 2.4, 0.55e-6, GRID_N; beam_diam_fraction=0.5f0)

    wf_e64 = cuda_wavefront_begin(Float64, 2.4, 0.55e-6, GRID_N; beam_diam_fraction=0.5)
    wf_e32 = cuda_wavefront_begin(Float32, 2.4, 0.55e-6, GRID_N; beam_diam_fraction=0.5f0)
    prop_circular_aperture(wf_e64, 0.6)
    prop_circular_aperture(wf_e32, 0.6f0)
    out_end64 = cuda_zeros(Float64, GRID_N, GRID_N)
    out_end32 = cuda_zeros(Float32, GRID_N, GRID_N)

    snap_q64 = capture_wavefront_state(wf_q64)
    snap_q32 = capture_wavefront_state(wf_q32)
    snap_p64 = capture_wavefront_state(wf_p64)
    snap_p32 = capture_wavefront_state(wf_p32)
    snap_w64 = capture_wavefront_state(wf_w64)
    snap_w32 = capture_wavefront_state(wf_w32)
    snap_s64 = capture_wavefront_state(wf_s64)
    snap_s32 = capture_wavefront_state(wf_s32)
    snap_a64 = capture_wavefront_state(wf_a64)
    snap_a32 = capture_wavefront_state(wf_a32)
    snap_e64 = capture_wavefront_state(wf_e64)
    snap_e32 = capture_wavefront_state(wf_e32)
    cuda_sync()

    # Warmup: exclude compilation and initial CUDA context setup.
    prop_qphase(wf_q64, 10.0, ctx_q64)
    prop_qphase(wf_q32, 10.0f0, ctx_q32)
    wf_p64.reference_surface = Proper.PLANAR
    prop_ptp(wf_p64, 0.01, ctx_p64)
    wf_p32.reference_surface = Proper.PLANAR
    prop_ptp(wf_p32, 0.01f0, ctx_p32)
    prop_wts(wf_w64, 0.01, ctx_w64)
    prop_wts(wf_w32, 0.01f0, ctx_w32)
    prop_stw(wf_s64, 0.01, ctx_s64)
    prop_stw(wf_s32, 0.01f0, ctx_s32)
    prop_circular_aperture(wf_a64, 0.6)
    prop_circular_aperture(wf_a32, 0.6f0)
    prop_end!(out_end64, wf_e64)
    prop_end!(out_end32, wf_e32)
    cuda_sync()

    q64 = run(@benchmarkable begin
        prop_qphase($wf_q64, 10.0, $ctx_q64)
        cuda_sync()
    end setup=(restore_wavefront_state!($wf_q64, $snap_q64); cuda_sync()) evals=1 samples=N_SAMPLES)

    q32 = run(@benchmarkable begin
        prop_qphase($wf_q32, 10.0f0, $ctx_q32)
        cuda_sync()
    end setup=(restore_wavefront_state!($wf_q32, $snap_q32); cuda_sync()) evals=1 samples=N_SAMPLES)

    p64 = run(@benchmarkable begin
        prop_ptp($wf_p64, 0.01, $ctx_p64)
        cuda_sync()
    end setup=(restore_wavefront_state!($wf_p64, $snap_p64); cuda_sync()) evals=1 samples=N_SAMPLES)

    p32 = run(@benchmarkable begin
        prop_ptp($wf_p32, 0.01f0, $ctx_p32)
        cuda_sync()
    end setup=(restore_wavefront_state!($wf_p32, $snap_p32); cuda_sync()) evals=1 samples=N_SAMPLES)

    w64 = run(@benchmarkable begin
        prop_wts($wf_w64, 0.01, $ctx_w64)
        cuda_sync()
    end setup=(restore_wavefront_state!($wf_w64, $snap_w64); cuda_sync()) evals=1 samples=N_SAMPLES)

    w32 = run(@benchmarkable begin
        prop_wts($wf_w32, 0.01f0, $ctx_w32)
        cuda_sync()
    end setup=(restore_wavefront_state!($wf_w32, $snap_w32); cuda_sync()) evals=1 samples=N_SAMPLES)

    s64 = run(@benchmarkable begin
        prop_stw($wf_s64, 0.01, $ctx_s64)
        cuda_sync()
    end setup=(restore_wavefront_state!($wf_s64, $snap_s64); cuda_sync()) evals=1 samples=N_SAMPLES)

    s32 = run(@benchmarkable begin
        prop_stw($wf_s32, 0.01f0, $ctx_s32)
        cuda_sync()
    end setup=(restore_wavefront_state!($wf_s32, $snap_s32); cuda_sync()) evals=1 samples=N_SAMPLES)

    a64 = run(@benchmarkable begin
        prop_circular_aperture($wf_a64, 0.6)
        cuda_sync()
    end setup=(restore_wavefront_state!($wf_a64, $snap_a64); cuda_sync()) evals=1 samples=N_SAMPLES)

    a32 = run(@benchmarkable begin
        prop_circular_aperture($wf_a32, 0.6f0)
        cuda_sync()
    end setup=(restore_wavefront_state!($wf_a32, $snap_a32); cuda_sync()) evals=1 samples=N_SAMPLES)

    e64 = run(@benchmarkable begin
        prop_end!($out_end64, $wf_e64)
        cuda_sync()
    end setup=(restore_wavefront_state!($wf_e64, $snap_e64); cuda_sync()) evals=1 samples=N_SAMPLES)

    e32 = run(@benchmarkable begin
        prop_end!($out_end32, $wf_e32)
        cuda_sync()
    end setup=(restore_wavefront_state!($wf_e32, $snap_e32); cuda_sync()) evals=1 samples=N_SAMPLES)

    report = Dict(
        "meta" => merge(cuda_report_meta(RUN_TAG; device=cuda_device_label()), Dict("grid_n" => GRID_N)),
        "policy" => "CUDA precision-split timing for steady-state propagation-heavy workload and selected kernels with per-sample wavefront state restore; TTFx excluded; per-sample synchronization included",
        "kernels" => Dict(
            "prop_qphase" => Dict("fp64" => trial_stats(q64), "fp32" => trial_stats(q32)),
            "prop_ptp" => Dict("fp64" => trial_stats(p64), "fp32" => trial_stats(p32)),
            "prop_wts" => Dict("fp64" => trial_stats(w64), "fp32" => trial_stats(w32)),
            "prop_stw" => Dict("fp64" => trial_stats(s64), "fp32" => trial_stats(s32)),
            "prop_circular_aperture" => Dict("fp64" => trial_stats(a64), "fp32" => trial_stats(a32)),
            "prop_end_mutating" => Dict("fp64" => trial_stats(e64), "fp32" => trial_stats(e32)),
        ),
    )

    return write_benchmark_report(REPORT_PATH, report)
end

benchmark_precision_split()
