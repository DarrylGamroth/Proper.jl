using BenchmarkTools
using JSON3
using Proper
include(joinpath(@__DIR__, "..", "..", "common", "metadata.jl"))
using .BenchMetadata

const RUN_TAG = "steady_state_supported_kernels"
const REPORT_PATH = joinpath(@__DIR__, "..", "..", "reports", "julia_supported_kernels.json")

@inline function trial_stats(t::BenchmarkTools.Trial)
    est = median(t)
    return Dict(
        "median_ns" => est.time,
        "median_allocs" => est.allocs,
        "median_bytes" => est.memory,
        "samples" => length(t.times),
    )
end

function bench_supported_kernels()
    nprop = 512
    nmap = 256
    samples = 20

    wf_q = prop_begin(2.4, 0.55e-6, nprop; beam_diam_fraction=0.5)
    ctx_q = RunContext(typeof(wf_q.field))

    wf_p = prop_begin(2.4, 0.55e-6, nprop; beam_diam_fraction=0.5)
    ctx_p = RunContext(typeof(wf_p.field))

    wf_a = prop_begin(2.4, 0.55e-6, nprop; beam_diam_fraction=0.5)
    out_end = zeros(Float64, nprop, nprop)

    img = rand(Float32, nmap, nmap)
    ctx_img = RunContext(typeof(img))
    rot_out = similar(img)
    mag_out = similar(img)
    szoom_out = similar(img)
    pix_out = similar(img, nmap ÷ 2, nmap ÷ 2)

    wf_map = prop_begin(1.0, 0.55e-6, nmap)
    dmap = rand(Float32, nmap, nmap)
    ctx_map = RunContext(typeof(dmap))
    res_out = zeros(Float64, nmap, nmap)
    res_opts = Proper.ResampleMapOptions(wf_map, wf_map.sampling_m, nmap / 2, nmap / 2)
    rect_out = zeros(Float64, nmap, nmap)
    round_out = similar(rect_out)

    # Warmup: exclude compilation from the reported timings.
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

    q = run(@benchmarkable begin
        prop_qphase($wf_q, 10.0, $ctx_q)
    end evals=1 samples=samples)

    p = run(@benchmarkable begin
        $wf_p.reference_surface = Proper.PLANAR
        prop_ptp($wf_p, 0.01, $ctx_p)
    end evals=1 samples=samples)

    a = run(@benchmarkable begin
        prop_circular_aperture($wf_a, 0.6)
    end evals=1 samples=samples)

    e = run(@benchmarkable begin
        prop_end!($out_end, $wf_a)
    end evals=1 samples=samples)

    r = run(@benchmarkable begin
        prop_rotate!($rot_out, $img, 12.0, $ctx_img)
    end evals=1 samples=samples)

    m = run(@benchmarkable begin
        prop_magnify!($mag_out, $img, 1.1, $ctx_img; QUICK=true)
    end evals=1 samples=samples)

    sz = run(@benchmarkable begin
        prop_szoom!($szoom_out, $img, 1.1)
    end evals=1 samples=samples)

    px = run(@benchmarkable begin
        prop_pixellate!($pix_out, $img, 2)
    end evals=1 samples=samples)

    rs = run(@benchmarkable begin
        prop_resamplemap!($res_out, $wf_map, $dmap, $res_opts, $ctx_map)
    end evals=1 samples=samples)

    rc = run(@benchmarkable begin
        prop_rectangle!($rect_out, $wf_map, 0.4, 0.2, 0.03, -0.05; ROTATION=22.0, NORM=true)
    end evals=1 samples=samples)

    rr = run(@benchmarkable begin
        prop_rounded_rectangle!($round_out, $wf_map, 0.05, 0.3, 0.2, 0.01, -0.02)
    end evals=1 samples=samples)

    return Dict(
        "meta" => merge(benchmark_metadata(run_tag=RUN_TAG), Dict("prop_grid_n" => nprop, "map_grid_n" => nmap)),
        "policy" => "steady-state supported CPU kernel timings; TTFx excluded",
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

report = bench_supported_kernels()
mkpath(dirname(REPORT_PATH))
open(REPORT_PATH, "w") do io
    JSON3.write(io, report)
end
println(report)
