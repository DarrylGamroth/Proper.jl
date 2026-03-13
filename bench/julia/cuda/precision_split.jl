using BenchmarkTools
using JSON3
using Proper
using CUDA
include(joinpath(@__DIR__, "..", "..", "common", "metadata.jl"))
include(joinpath(@__DIR__, "..", "..", "common", "cuda_support.jl"))
using .BenchMetadata

const RUN_TAG = "steady_state_cuda_precision_split"
const REPORT_PATH = joinpath(@__DIR__, "..", "..", "reports", "cuda_precision_split.json")
const GRID_N = 512
const N_SAMPLES = 20

CUDA.allowscalar(false)

function steady_state_workload(::Type{T}) where {T<:AbstractFloat}
    wf = cuda_wavefront_begin(T, 2.4, 0.55e-6, GRID_N; beam_diam_fraction=T(0.5))
    ctx = RunContext(wf)
    prop_circular_aperture(wf, T(0.6))
    prop_lens(wf, T(20.0))
    prop_propagate(wf, T(20.0), ctx)
    prop_end(wf)
    cuda_sync()
    return nothing
end

function benchmark_precision_split()
    wf_q64 = cuda_wavefront_begin(Float64, 2.4, 0.55e-6, GRID_N; beam_diam_fraction=0.5)
    ctx_q64 = RunContext(wf_q64)
    wf_q32 = cuda_wavefront_begin(Float32, 2.4, 0.55e-6, GRID_N; beam_diam_fraction=0.5f0)
    ctx_q32 = RunContext(wf_q32)

    wf_p64 = cuda_wavefront_begin(Float64, 2.4, 0.55e-6, GRID_N; beam_diam_fraction=0.5)
    ctx_p64 = RunContext(wf_p64)
    wf_p32 = cuda_wavefront_begin(Float32, 2.4, 0.55e-6, GRID_N; beam_diam_fraction=0.5f0)
    ctx_p32 = RunContext(wf_p32)

    wf_a64 = cuda_wavefront_begin(Float64, 2.4, 0.55e-6, GRID_N; beam_diam_fraction=0.5)
    wf_a32 = cuda_wavefront_begin(Float32, 2.4, 0.55e-6, GRID_N; beam_diam_fraction=0.5f0)

    out_end64 = cuda_zeros(Float64, GRID_N, GRID_N)
    out_end32 = cuda_zeros(Float32, GRID_N, GRID_N)

    # Warmup: exclude compilation and initial CUDA context setup.
    steady_state_workload(Float64)
    steady_state_workload(Float32)
    prop_qphase(wf_q64, 10.0, ctx_q64)
    prop_qphase(wf_q32, 10.0f0, ctx_q32)
    wf_p64.reference_surface = Proper.PLANAR
    prop_ptp(wf_p64, 0.01, ctx_p64)
    wf_p32.reference_surface = Proper.PLANAR
    prop_ptp(wf_p32, 0.01f0, ctx_p32)
    prop_circular_aperture(wf_a64, 0.6)
    prop_circular_aperture(wf_a32, 0.6f0)
    prop_end!(out_end64, wf_a64)
    prop_end!(out_end32, wf_a32)
    cuda_sync()

    workload64 = run(@benchmarkable begin
        steady_state_workload(Float64)
        cuda_sync()
    end evals=1 samples=N_SAMPLES)

    workload32 = run(@benchmarkable begin
        steady_state_workload(Float32)
        cuda_sync()
    end evals=1 samples=N_SAMPLES)

    q64 = run(@benchmarkable begin
        prop_qphase($wf_q64, 10.0, $ctx_q64)
        cuda_sync()
    end evals=1 samples=N_SAMPLES)

    q32 = run(@benchmarkable begin
        prop_qphase($wf_q32, 10.0f0, $ctx_q32)
        cuda_sync()
    end evals=1 samples=N_SAMPLES)

    p64 = run(@benchmarkable begin
        $wf_p64.reference_surface = Proper.PLANAR
        prop_ptp($wf_p64, 0.01, $ctx_p64)
        cuda_sync()
    end evals=1 samples=N_SAMPLES)

    p32 = run(@benchmarkable begin
        $wf_p32.reference_surface = Proper.PLANAR
        prop_ptp($wf_p32, 0.01f0, $ctx_p32)
        cuda_sync()
    end evals=1 samples=N_SAMPLES)

    a64 = run(@benchmarkable begin
        prop_circular_aperture($wf_a64, 0.6)
        cuda_sync()
    end evals=1 samples=N_SAMPLES)

    a32 = run(@benchmarkable begin
        prop_circular_aperture($wf_a32, 0.6f0)
        cuda_sync()
    end evals=1 samples=N_SAMPLES)

    e64 = run(@benchmarkable begin
        prop_end!($out_end64, $wf_a64)
        cuda_sync()
    end evals=1 samples=N_SAMPLES)

    e32 = run(@benchmarkable begin
        prop_end!($out_end32, $wf_a32)
        cuda_sync()
    end evals=1 samples=N_SAMPLES)

    report = Dict(
        "meta" => merge(cuda_report_meta(RUN_TAG; device=cuda_device_label()), Dict("grid_n" => GRID_N)),
        "policy" => "CUDA precision-split timing for steady-state propagation-heavy workload and selected kernels; TTFx excluded; per-sample synchronization included",
        "workloads" => Dict(
            "steady_state_fp64" => trial_stats(workload64),
            "steady_state_fp32" => trial_stats(workload32),
        ),
        "kernels" => Dict(
            "prop_qphase" => Dict("fp64" => trial_stats(q64), "fp32" => trial_stats(q32)),
            "prop_ptp" => Dict("fp64" => trial_stats(p64), "fp32" => trial_stats(p32)),
            "prop_circular_aperture" => Dict("fp64" => trial_stats(a64), "fp32" => trial_stats(a32)),
            "prop_end_mutating" => Dict("fp64" => trial_stats(e64), "fp32" => trial_stats(e32)),
        ),
    )

    return write_benchmark_report(REPORT_PATH, report)
end

benchmark_precision_split()
