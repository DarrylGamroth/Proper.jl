using BenchmarkTools
using JSON3
using Proper

function bench_phase2_kernels()
    wf_q = prop_begin(2.4, 0.55e-6, 512; beam_diam_fraction=0.5)
    wf_l = prop_begin(2.4, 0.55e-6, 512; beam_diam_fraction=0.5)
    wf_p = prop_begin(2.4, 0.55e-6, 512; beam_diam_fraction=0.5)

    # Warmup (exclude compilation)
    prop_qphase(wf_q, 10.0)
    prop_lens(wf_l, 20.0)
    wf_p.reference_surface = :PLANAR
    prop_ptp(wf_p, 0.01)

    q = @benchmarkable prop_qphase($wf_q, 10.0)
    l = @benchmarkable prop_lens($wf_l, 20.0)
    p = @benchmarkable begin
        $wf_p.reference_surface = :PLANAR
        prop_ptp($wf_p, 0.01)
    end

    qr = run(q)
    lr = run(l)
    pr = run(p)

    return Dict(
        "prop_qphase_median_ns" => median(qr).time,
        "prop_lens_median_ns" => median(lr).time,
        "prop_ptp_median_ns" => median(pr).time,
    )
end

report = bench_phase2_kernels()
out = joinpath(@__DIR__, "..", "..", "reports", "phase2_kernels.json")
mkpath(dirname(out))
open(out, "w") do io
    JSON3.write(io, report)
end
println(report)
