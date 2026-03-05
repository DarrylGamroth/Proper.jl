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

function bench_phase2_kernels()
    wf_q = prop_begin(2.4, 0.55e-6, 512; beam_diam_fraction=0.5)
    wf_l = prop_begin(2.4, 0.55e-6, 512; beam_diam_fraction=0.5)
    wf_p = prop_begin(2.4, 0.55e-6, 512; beam_diam_fraction=0.5)
    ctx_p = RunContext(typeof(wf_p.field))

    # Warmup (exclude compilation)
    prop_qphase(wf_q, 10.0)
    prop_lens(wf_l, 20.0)
    wf_p.reference_surface = Proper.PLANAR
    prop_ptp(wf_p, 0.01, ctx_p)

    q = @benchmarkable prop_qphase($wf_q, 10.0) evals=1 samples=40
    l = @benchmarkable prop_lens($wf_l, 20.0) evals=1 samples=40
    p = @benchmarkable begin
        $wf_p.reference_surface = Proper.PLANAR
        prop_ptp($wf_p, 0.01, $ctx_p)
    end evals=1 samples=40

    qr = trial_stats(run(q))
    lr = trial_stats(run(l))
    pr = trial_stats(run(p))

    return Dict(
        "meta" => benchmark_metadata(run_tag="steady_state_kernels_phase2"),
        "policy" => "steady-state kernel timing only; TTFx excluded",
        "kernels" => Dict(
            "prop_qphase" => qr,
            "prop_lens" => lr,
            "prop_ptp" => pr,
        ),
    )
end

report = bench_phase2_kernels()
out = joinpath(@__DIR__, "..", "..", "reports", "phase2_kernels.json")
mkpath(dirname(out))
open(out, "w") do io
    JSON3.write(io, report)
end
println(report)
