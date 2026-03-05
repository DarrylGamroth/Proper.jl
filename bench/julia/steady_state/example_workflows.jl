using BenchmarkTools
using JSON3
using Proper
include(joinpath(@__DIR__, "..", "..", "common", "metadata.jl"))
using .BenchMetadata

const EXAMPLES_DIR = normpath(joinpath(@__DIR__, "..", "..", "..", "examples"))
include(joinpath(EXAMPLES_DIR, "simple_prescription.jl"))
include(joinpath(EXAMPLES_DIR, "simple_telescope.jl"))
include(joinpath(EXAMPLES_DIR, "psdtest.jl"))

@inline function trial_stats(t::BenchmarkTools.Trial)
    est = median(t)
    return Dict(
        "median_ns" => est.time,
        "median_allocs" => est.allocs,
        "median_bytes" => est.memory,
        "samples" => length(t.times),
    )
end

function bench_example_workflows()
    λ = 0.55e-6

    # Warmup (exclude compilation)
    simple_prescription(λ, 256)
    simple_telescope(λ, 256)
    psdtest(λ, 128, Dict("usepsdmap" => true))

    sp = run(@benchmarkable simple_prescription($λ, 256) evals=1 samples=20)
    st = run(@benchmarkable simple_telescope($λ, 256) evals=1 samples=20)
    ps = run(@benchmarkable psdtest($λ, 128, Dict("usepsdmap" => true)) evals=1 samples=10)

    return Dict(
        "meta" => benchmark_metadata(run_tag="steady_state_examples"),
        "policy" => "steady-state example timing only; TTFx excluded",
        "examples" => Dict(
            "simple_prescription_256" => trial_stats(sp),
            "simple_telescope_256" => trial_stats(st),
            "psdtest_128" => trial_stats(ps),
        ),
    )
end

report = bench_example_workflows()
out = joinpath(@__DIR__, "..", "..", "reports", "example_workflows.json")
mkpath(dirname(out))
open(out, "w") do io
    JSON3.write(io, report)
end
println(report)
