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

function bench_ka_interp_kernels()
    n = 256
    a = reshape(collect(Float32, 1:(n * n)), n, n)
    out_loop = similar(a)
    out_ka = similar(a)
    x = collect(Float32, 1:n)
    y = collect(Float32, 1:n)
    sty = Proper.CubicInterpStyle()

    img = rand(Float32, n, n)
    rot_loop = similar(img)
    rot_ka = similar(img)
    ang = -deg2rad(12.0)
    c = cos(ang)
    s = sin(ang)
    opts_cubic = Proper.RotateOptions(img, pairs((; CUBIC=true)))
    opts_linear = Proper.RotateOptions(img, pairs((; METH="linear")))

    wf = prop_begin(1.0, 0.55e-6, n)
    ctx = RunContext(typeof(wf.field))
    dmap = rand(Float32, n, n)
    res_out = zeros(Float64, n, n)
    res_opts = Proper.ResampleMapOptions(wf, wf.sampling_m, n / 2, n / 2)
    mag_out = similar(img)

    Proper._prop_cubic_conv_grid_loop!(out_loop, sty, a, x, y)
    Proper.ka_cubic_conv_grid!(out_ka, a, x, y)
    Proper._prop_rotate_cubic!(sty, rot_loop, img, c, s, opts_cubic)
    Proper.ka_rotate_cubic!(rot_ka, img, c, s, opts_cubic.cx, opts_cubic.cy, opts_cubic.sx, opts_cubic.sy)
    Proper._prop_rotate_linear!(rot_loop, img, c, s, opts_linear)
    Proper.ka_rotate_linear!(rot_ka, img, c, s, opts_linear.cx, opts_linear.cy, opts_linear.sx, opts_linear.sy)
    prop_resamplemap!(res_out, wf, dmap, res_opts, ctx)
    prop_magnify!(mag_out, img, 1.1, ctx; QUICK=true)

    samples = 20
    cubic_pair = pair_stats(
        run(@benchmarkable Proper._prop_cubic_conv_grid_loop!($out_loop, $sty, $a, $x, $y) evals=1 samples=samples),
        run(@benchmarkable Proper.ka_cubic_conv_grid!($out_ka, $a, $x, $y) evals=1 samples=samples),
    )
    rotate_cubic_pair = pair_stats(
        run(@benchmarkable Proper._prop_rotate_cubic!($sty, $rot_loop, $img, $c, $s, $opts_cubic) evals=1 samples=samples),
        run(@benchmarkable Proper.ka_rotate_cubic!($rot_ka, $img, $c, $s, $(opts_cubic.cx), $(opts_cubic.cy), $(opts_cubic.sx), $(opts_cubic.sy)) evals=1 samples=samples),
    )
    rotate_linear_pair = pair_stats(
        run(@benchmarkable Proper._prop_rotate_linear!($rot_loop, $img, $c, $s, $opts_linear) evals=1 samples=samples),
        run(@benchmarkable Proper.ka_rotate_linear!($rot_ka, $img, $c, $s, $(opts_linear.cx), $(opts_linear.cy), $(opts_linear.sx), $(opts_linear.sy)) evals=1 samples=samples),
    )

    resample_trial = run(@benchmarkable prop_resamplemap!($res_out, $wf, $dmap, $res_opts, $ctx) evals=1 samples=samples)
    magnify_trial = run(@benchmarkable prop_magnify!($mag_out, $img, 1.1, $ctx; QUICK=true) evals=1 samples=samples)

    return Dict(
        "meta" => merge(benchmark_metadata(run_tag="steady_state_ka_interp"), Dict("grid_n" => n)),
        "policy" => "KA interpolation pilot benchmark; loop vs KA internal kernels plus current public wrapper timings",
        "pairs" => Dict(
            "cubic_conv_grid" => cubic_pair,
            "rotate_cubic" => rotate_cubic_pair,
            "rotate_linear" => rotate_linear_pair,
        ),
        "public" => Dict(
            "prop_resamplemap_mutating" => trial_stats(resample_trial),
            "prop_magnify_quick_mutating" => trial_stats(magnify_trial),
        ),
    )
end

report = bench_ka_interp_kernels()
out = joinpath(@__DIR__, "..", "..", "reports", "ka_interp_kernels.json")
mkpath(dirname(out))
open(out, "w") do io
    JSON3.write(io, report)
end
println(report)
