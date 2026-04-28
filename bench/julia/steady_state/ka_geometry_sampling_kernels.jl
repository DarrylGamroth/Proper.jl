using BenchmarkTools
using JSON3
using Proper
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

@inline function pair_stats(loop_trial::BenchmarkTools.Trial, ka_trial::BenchmarkTools.Trial)
    loop_stats = trial_stats(loop_trial)
    ka_stats = trial_stats(ka_trial)
    loop_ns = Float64(loop_stats["median_ns"])
    ka_ns = Float64(ka_stats["median_ns"])
    return Dict(
        "loop" => loop_stats,
        "ka" => ka_stats,
        "speedup_loop_over_ka" => loop_ns / max(ka_ns, 1.0),
    )
end

function bench_ka_geometry_sampling_kernels()
    samples = 20
    n = 256
    wf = prop_begin(1.0, 0.55e-6, n)
    RT = real(eltype(wf.field))

    ellipse_loop = zeros(RT, n, n)
    ellipse_ka = similar(ellipse_loop)
    rect_loop = zeros(RT, n, n)
    rect_ka = similar(rect_loop)
    xverts = RT[-0.20, 0.12, 0.28, -0.08]
    yverts = RT[-0.18, -0.22, 0.19, 0.25]
    ipoly_loop = zeros(RT, n, n)
    ipoly_ka = similar(ipoly_loop)
    round_loop = zeros(RT, n, n)
    round_ka = similar(round_loop)

    img = rand(Float32, 128, 128)
    pix_loop = similar(img, 64, 64)
    pix_ka = similar(pix_loop)
    szoom_loop = similar(img, 96, 96)
    szoom_ka = similar(szoom_loop)
    szoom_loop_ws = Proper.SamplingWorkspace(typeof(szoom_loop), Float32)
    szoom_ka_ws = Proper.SamplingWorkspace(typeof(szoom_ka), Float32)

    Proper._prop_ellipse!(Proper.GeometryLoopExecStyle(), ellipse_loop, wf, 0.35, 0.25, 0.05, -0.03, Proper.EllipseOptions(pairs((; ROTATION=13.0, DARK=true, NORM=true))))
    Proper._prop_ellipse!(Proper.GeometryKAExecStyle(), ellipse_ka, wf, 0.35, 0.25, 0.05, -0.03, Proper.EllipseOptions(pairs((; ROTATION=13.0, DARK=true, NORM=true))))
    Proper._prop_rectangle!(Proper.GeometryLoopExecStyle(), rect_loop, wf, 0.4, 0.2, 0.03, -0.05, Proper.RectangleOptions(pairs((; ROTATION=22.0, NORM=true))))
    Proper._prop_rectangle!(Proper.GeometryKAExecStyle(), rect_ka, wf, 0.4, 0.2, 0.03, -0.05, Proper.RectangleOptions(pairs((; ROTATION=22.0, NORM=true))))
    Proper._prop_irregular_polygon(Proper.GeometryLoopExecStyle(), ipoly_loop, wf, xverts, yverts, Proper.IrregularPolygonOptions(pairs((; NORM=true))))
    Proper._prop_irregular_polygon(Proper.GeometryKAExecStyle(), ipoly_ka, wf, xverts, yverts, Proper.IrregularPolygonOptions(pairs((; NORM=true))))
    Proper._prop_rounded_rectangle!(Proper.GeometryLoopExecStyle(), round_loop, wf, 0.05, 0.3, 0.2, 0.01, -0.02)
    Proper._prop_rounded_rectangle!(Proper.GeometryKAExecStyle(), round_ka, wf, 0.05, 0.3, 0.2, 0.01, -0.02)
    Proper._prop_pixellate_factor!(Proper.SamplingLoopExecStyle(), pix_loop, img, 2)
    Proper._prop_pixellate_factor!(Proper.SamplingKAExecStyle(), pix_ka, img, 2)
    Proper._prop_szoom!(Proper.SamplingLoopExecStyle(), szoom_loop, img, Float32(1.35), szoom_loop_ws)
    Proper._prop_szoom!(Proper.SamplingKAExecStyle(), szoom_ka, img, Float32(1.35), szoom_ka_ws)

    ellipse_pair = pair_stats(
        run(@benchmarkable Proper._prop_ellipse!(Proper.GeometryLoopExecStyle(), $ellipse_loop, $wf, 0.35, 0.25, 0.05, -0.03, $(Proper.EllipseOptions(pairs((; ROTATION=13.0, DARK=true, NORM=true))))) evals=1 samples=samples),
        run(@benchmarkable Proper._prop_ellipse!(Proper.GeometryKAExecStyle(), $ellipse_ka, $wf, 0.35, 0.25, 0.05, -0.03, $(Proper.EllipseOptions(pairs((; ROTATION=13.0, DARK=true, NORM=true))))) evals=1 samples=samples),
    )
    rectangle_pair = pair_stats(
        run(@benchmarkable Proper._prop_rectangle!(Proper.GeometryLoopExecStyle(), $rect_loop, $wf, 0.4, 0.2, 0.03, -0.05, $(Proper.RectangleOptions(pairs((; ROTATION=22.0, NORM=true))))) evals=1 samples=samples),
        run(@benchmarkable Proper._prop_rectangle!(Proper.GeometryKAExecStyle(), $rect_ka, $wf, 0.4, 0.2, 0.03, -0.05, $(Proper.RectangleOptions(pairs((; ROTATION=22.0, NORM=true))))) evals=1 samples=samples),
    )
    ipoly_pair = pair_stats(
        run(@benchmarkable Proper._prop_irregular_polygon(Proper.GeometryLoopExecStyle(), $ipoly_loop, $wf, $xverts, $yverts, $(Proper.IrregularPolygonOptions(pairs((; NORM=true))))) evals=1 samples=samples),
        run(@benchmarkable Proper._prop_irregular_polygon(Proper.GeometryKAExecStyle(), $ipoly_ka, $wf, $xverts, $yverts, $(Proper.IrregularPolygonOptions(pairs((; NORM=true))))) evals=1 samples=samples),
    )
    round_pair = pair_stats(
        run(@benchmarkable Proper._prop_rounded_rectangle!(Proper.GeometryLoopExecStyle(), $round_loop, $wf, 0.05, 0.3, 0.2, 0.01, -0.02) evals=1 samples=samples),
        run(@benchmarkable Proper._prop_rounded_rectangle!(Proper.GeometryKAExecStyle(), $round_ka, $wf, 0.05, 0.3, 0.2, 0.01, -0.02) evals=1 samples=samples),
    )
    pix_pair = pair_stats(
        run(@benchmarkable Proper._prop_pixellate_factor!(Proper.SamplingLoopExecStyle(), $pix_loop, $img, 2) evals=1 samples=samples),
        run(@benchmarkable Proper._prop_pixellate_factor!(Proper.SamplingKAExecStyle(), $pix_ka, $img, 2) evals=1 samples=samples),
    )
    szoom_pair = pair_stats(
        run(@benchmarkable Proper._prop_szoom!(Proper.SamplingLoopExecStyle(), $szoom_loop, $img, Float32(1.35), $szoom_loop_ws) evals=1 samples=samples),
        run(@benchmarkable Proper._prop_szoom!(Proper.SamplingKAExecStyle(), $szoom_ka, $img, Float32(1.35), $szoom_ka_ws) evals=1 samples=samples),
    )

    return Dict(
        "meta" => merge(benchmark_metadata(run_tag="steady_state_ka_geometry_sampling"), Dict("grid_n" => n, "sample_grid_n" => 128)),
        "policy" => "KA geometry/sampling pilot benchmark; loop vs KA internal kernels",
        "pairs" => Dict(
            "ellipse" => ellipse_pair,
            "rectangle" => rectangle_pair,
            "irregular_polygon" => ipoly_pair,
            "rounded_rectangle" => round_pair,
            "pixellate" => pix_pair,
            "szoom" => szoom_pair,
        ),
    )
end

report = bench_ka_geometry_sampling_kernels()
out = joinpath(@__DIR__, "..", "..", "reports", "ka_geometry_sampling_kernels.json")
mkpath(dirname(out))
open(out, "w") do io
    JSON3.write(io, report)
end
println(report)
