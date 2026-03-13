using BenchmarkTools
using JSON3
using Proper
include(joinpath(@__DIR__, "..", "..", "common", "metadata.jl"))
include(joinpath(@__DIR__, "..", "..", "common", "cuda_support.jl"))
using .BenchMetadata

const RUN_TAG = "steady_state_cuda_supported_kernels"
const REPORT_PATH = joinpath(@__DIR__, "..", "..", "reports", "cuda_supported_kernels.json")

function bench_cuda_supported_kernels(cuda_mod)
    nprop = 512
    nmap = 256
    samples = 20

    wf_q = cuda_wavefront_begin(cuda_mod, 2.4, 0.55e-6, nprop; beam_diam_fraction=0.5)
    ctx_q = RunContext(typeof(wf_q.field))

    wf_p = cuda_wavefront_begin(cuda_mod, 2.4, 0.55e-6, nprop; beam_diam_fraction=0.5)
    ctx_p = RunContext(typeof(wf_p.field))

    wf_a = cuda_wavefront_begin(cuda_mod, 2.4, 0.55e-6, nprop; beam_diam_fraction=0.5)
    out_end = cuda_mod.zeros(Float64, nprop, nprop)

    img = cuda_mod.rand(Float32, nmap, nmap)
    ctx_img = RunContext(typeof(img))
    rot_out = similar(img)
    mag_out = similar(img)
    szoom_out = similar(img)
    pix_out = similar(img, nmap ÷ 2, nmap ÷ 2)

    wf_map = cuda_wavefront_begin(cuda_mod, 1.0, 0.55e-6, nmap)
    dmap = cuda_mod.rand(Float32, nmap, nmap)
    ctx_map = RunContext(typeof(dmap))
    res_out = cuda_mod.zeros(Float64, nmap, nmap)
    res_opts = Proper.ResampleMapOptions(wf_map, wf_map.sampling_m, nmap / 2, nmap / 2)
    rect_out = cuda_mod.zeros(Float64, nmap, nmap)
    round_out = similar(rect_out)

    # Warmup
    prop_qphase(wf_q, 10.0, ctx_q)
    wf_p.reference_surface = Proper.PLANAR
    prop_ptp(wf_p, 0.01, ctx_p)
    prop_circular_aperture(wf_a, 0.6)
    prop_end!(out_end, wf_a)
    prop_rotate!(rot_out, img, 12.0, ctx_img)
    prop_magnify!(mag_out, img, 1.1, ctx_img; QUICK=true)
    prop_szoom!(szoom_out, img, 1.1)
    prop_pixellate!(pix_out, img, 2)
    prop_resamplemap!(res_out, wf_map, dmap, res_opts, ctx_map)
    prop_rectangle!(rect_out, wf_map, 0.4, 0.2, 0.03, -0.05; ROTATION=22.0, NORM=true)
    prop_rounded_rectangle!(round_out, wf_map, 0.05, 0.3, 0.2, 0.01, -0.02)
    cuda_mod.synchronize()

    q = run(@benchmarkable begin
        prop_qphase($wf_q, 10.0, $ctx_q)
        $cuda_mod.synchronize()
    end evals=1 samples=samples)

    p = run(@benchmarkable begin
        $wf_p.reference_surface = Proper.PLANAR
        prop_ptp($wf_p, 0.01, $ctx_p)
        $cuda_mod.synchronize()
    end evals=1 samples=samples)

    a = run(@benchmarkable begin
        prop_circular_aperture($wf_a, 0.6)
        $cuda_mod.synchronize()
    end evals=1 samples=samples)

    e = run(@benchmarkable begin
        prop_end!($out_end, $wf_a)
        $cuda_mod.synchronize()
    end evals=1 samples=samples)

    r = run(@benchmarkable begin
        prop_rotate!($rot_out, $img, 12.0, $ctx_img)
        $cuda_mod.synchronize()
    end evals=1 samples=samples)

    m = run(@benchmarkable begin
        prop_magnify!($mag_out, $img, 1.1, $ctx_img; QUICK=true)
        $cuda_mod.synchronize()
    end evals=1 samples=samples)

    sz = run(@benchmarkable begin
        prop_szoom!($szoom_out, $img, 1.1)
        $cuda_mod.synchronize()
    end evals=1 samples=samples)

    px = run(@benchmarkable begin
        prop_pixellate!($pix_out, $img, 2)
        $cuda_mod.synchronize()
    end evals=1 samples=samples)

    rs = run(@benchmarkable begin
        prop_resamplemap!($res_out, $wf_map, $dmap, $res_opts, $ctx_map)
        $cuda_mod.synchronize()
    end evals=1 samples=samples)

    rc = run(@benchmarkable begin
        prop_rectangle!($rect_out, $wf_map, 0.4, 0.2, 0.03, -0.05; ROTATION=22.0, NORM=true)
        $cuda_mod.synchronize()
    end evals=1 samples=samples)

    rr = run(@benchmarkable begin
        prop_rounded_rectangle!($round_out, $wf_map, 0.05, 0.3, 0.2, 0.01, -0.02)
        $cuda_mod.synchronize()
    end evals=1 samples=samples)

    return Dict(
        "meta" => merge(cuda_report_meta(RUN_TAG, cuda_mod), Dict("prop_grid_n" => nprop, "map_grid_n" => nmap)),
        "policy" => "steady-state supported CUDA kernel timings; TTFx excluded; per-sample synchronization included",
        "kernels" => Dict(
            "prop_qphase" => trial_stats(q),
            "prop_ptp" => trial_stats(p),
            "prop_circular_aperture" => trial_stats(a),
            "prop_end_mutating" => trial_stats(e),
            "prop_rotate_mutating" => trial_stats(r),
            "prop_magnify_quick_mutating" => trial_stats(m),
            "prop_szoom_mutating" => trial_stats(sz),
            "prop_pixellate_mutating" => trial_stats(px),
            "prop_resamplemap_mutating" => trial_stats(rs),
            "prop_rectangle_mutating" => trial_stats(rc),
            "prop_rounded_rectangle_mutating" => trial_stats(rr),
        ),
    )
end

function main()
    cuda_mod, reason = load_cuda_backend()
    if cuda_mod === nothing
        return write_benchmark_report(REPORT_PATH, skipped_cuda_report(RUN_TAG, reason))
    end
    return write_benchmark_report(REPORT_PATH, bench_cuda_supported_kernels(cuda_mod))
end

main()
