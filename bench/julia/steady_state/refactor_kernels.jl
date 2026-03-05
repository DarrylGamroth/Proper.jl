using BenchmarkTools
using JSON3
using Proper
using Random
include(joinpath(@__DIR__, "..", "..", "common", "metadata.jl"))
using .BenchMetadata

@inline function trial_stats(t::BenchmarkTools.Trial)
    est = median(t)
    return Dict(
        "median_ns" => est.time,
        "median_allocs" => est.allocs,
        "median_bytes" => est.memory,
        "samples" => length(t.times),
    )
end

@inline function pair_stats(wrapper::BenchmarkTools.Trial, mutating::BenchmarkTools.Trial)
    w = trial_stats(wrapper)
    m = trial_stats(mutating)
    wn = Float64(w["median_ns"])
    mn = Float64(m["median_ns"])
    wb = Float64(w["median_bytes"])
    mb = Float64(m["median_bytes"])
    return Dict(
        "wrapper" => w,
        "mutating" => m,
        "speedup_wrapper_over_mutating" => wn / max(mn, 1.0),
        "median_byte_reduction" => wb - mb,
    )
end

function bench_refactor_kernels()
    rng = MersenneTwister(11)

    # Interpolation / map kernels
    wf_map = prop_begin(1.0, 0.55e-6, 128)
    ctx_map = RunContext(typeof(wf_map.field))
    dmap = rand(rng, Float32, 128, 128)
    res_out = zeros(Float64, 128, 128)
    res_opts = Proper.ResampleMapOptions(wf_map, wf_map.sampling_m, 64.0, 64.0)

    img = rand(rng, Float32, 128, 128)
    rot_out = similar(img)
    mag_out = similar(img, 128, 128)

    # Geometry kernels
    wf_geom = prop_begin(1.0, 0.55e-6, 128)
    RT = real(eltype(wf_geom.field))
    ellipse_out = zeros(RT, 128, 128)
    rect_out = zeros(RT, 128, 128)
    poly_out = zeros(RT, 128, 128)
    ipoly_out = zeros(RT, 128, 128)
    xverts = RT[-0.2, 0.12, 0.28, -0.08]
    yverts = RT[-0.18, -0.22, 0.19, 0.25]

    # PSD hotspot
    wf_psd = prop_begin(0.212, 0.55e-6, 128)

    # Warmup (exclude compilation)
    prop_resamplemap(wf_map, dmap, wf_map.sampling_m, 64.0, 64.0, ctx_map)
    prop_resamplemap!(res_out, wf_map, dmap, res_opts, ctx_map)
    prop_rotate(img, 12.0, ctx_map)
    prop_rotate!(rot_out, img, 12.0, ctx_map)
    prop_magnify(img, 1.1, 128, ctx_map; QUICK=true)
    prop_magnify!(mag_out, img, 1.1, ctx_map; QUICK=true)

    prop_ellipse(wf_geom, 0.35, 0.2, 0.0, 0.0; NORM=true)
    prop_ellipse!(ellipse_out, wf_geom, 0.35, 0.2, 0.0, 0.0; NORM=true)
    prop_rectangle(wf_geom, 0.4, 0.2, 0.0, 0.0; NORM=true)
    prop_rectangle!(rect_out, wf_geom, 0.4, 0.2, 0.0, 0.0; NORM=true)
    prop_polygon(wf_geom, 6, 0.33, 0.0, 0.0; NORM=true)
    prop_polygon!(poly_out, wf_geom, 6, 0.33, 0.0, 0.0; NORM=true)
    prop_irregular_polygon(wf_geom, xverts, yverts; NORM=true)
    prop_irregular_polygon!(ipoly_out, wf_geom, xverts, yverts; NORM=true)

    prop_psd_errormap(wf_psd, 3.29e-23, 212.26, 7.8; no_apply=true, rng=MersenneTwister(3))

    samples_fast = 30
    samples_heavy = 15

    resample_pair = pair_stats(
        run(@benchmarkable prop_resamplemap($wf_map, $dmap, $wf_map.sampling_m, 64.0, 64.0, $ctx_map) evals=1 samples=samples_fast),
        run(@benchmarkable prop_resamplemap!($res_out, $wf_map, $dmap, $res_opts, $ctx_map) evals=1 samples=samples_fast),
    )

    rotate_pair = pair_stats(
        run(@benchmarkable prop_rotate($img, 12.0, $ctx_map) evals=1 samples=samples_fast),
        run(@benchmarkable prop_rotate!($rot_out, $img, 12.0, $ctx_map) evals=1 samples=samples_fast),
    )

    magnify_pair = pair_stats(
        run(@benchmarkable prop_magnify($img, 1.1, 128, $ctx_map; QUICK=true) evals=1 samples=samples_fast),
        run(@benchmarkable prop_magnify!($mag_out, $img, 1.1, $ctx_map; QUICK=true) evals=1 samples=samples_fast),
    )

    ellipse_pair = pair_stats(
        run(@benchmarkable prop_ellipse($wf_geom, 0.35, 0.2, 0.0, 0.0; NORM=true) evals=1 samples=samples_fast),
        run(@benchmarkable prop_ellipse!($ellipse_out, $wf_geom, 0.35, 0.2, 0.0, 0.0; NORM=true) evals=1 samples=samples_fast),
    )

    rectangle_pair = pair_stats(
        run(@benchmarkable prop_rectangle($wf_geom, 0.4, 0.2, 0.0, 0.0; NORM=true) evals=1 samples=samples_fast),
        run(@benchmarkable prop_rectangle!($rect_out, $wf_geom, 0.4, 0.2, 0.0, 0.0; NORM=true) evals=1 samples=samples_fast),
    )

    polygon_pair = pair_stats(
        run(@benchmarkable prop_polygon($wf_geom, 6, 0.33, 0.0, 0.0; NORM=true) evals=1 samples=samples_fast),
        run(@benchmarkable prop_polygon!($poly_out, $wf_geom, 6, 0.33, 0.0, 0.0; NORM=true) evals=1 samples=samples_fast),
    )

    irregular_polygon_pair = pair_stats(
        run(@benchmarkable prop_irregular_polygon($wf_geom, $xverts, $yverts; NORM=true) evals=1 samples=samples_fast),
        run(@benchmarkable prop_irregular_polygon!($ipoly_out, $wf_geom, $xverts, $yverts; NORM=true) evals=1 samples=samples_fast),
    )

    psd_trial = run(@benchmarkable prop_psd_errormap($wf_psd, 3.29e-23, 212.26, 7.8; no_apply=true, rng=MersenneTwister(3)) evals=1 samples=samples_heavy)

    return Dict(
        "meta" => benchmark_metadata(run_tag="steady_state_kernels_refactor"),
        "policy" => "steady-state kernel timing only; TTFx excluded",
        "pairs" => Dict(
            "resamplemap" => resample_pair,
            "rotate" => rotate_pair,
            "magnify_quick" => magnify_pair,
            "ellipse" => ellipse_pair,
            "rectangle" => rectangle_pair,
            "polygon" => polygon_pair,
            "irregular_polygon" => irregular_polygon_pair,
        ),
        "hotspots" => Dict(
            "psd_errormap_no_apply" => trial_stats(psd_trial),
        ),
    )
end

report = bench_refactor_kernels()
out = joinpath(@__DIR__, "..", "..", "reports", "refactor_kernels.json")
mkpath(dirname(out))
open(out, "w") do io
    JSON3.write(io, report)
end
println(report)
