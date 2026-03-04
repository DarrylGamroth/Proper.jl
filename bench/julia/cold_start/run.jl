using JSON3
using proper
include(joinpath(@__DIR__, "..", "..", "common", "metadata.jl"))
using .BenchMetadata

function workload()
    wf = prop_begin(2.4, 0.55e-6, 512; beam_diam_fraction=0.5)
    prop_circular_aperture(wf, 0.6)
    prop_lens(wf, 20.0)
    prop_propagate(wf, 20.0)
    prop_end(wf)
end

start_ns = time_ns()
workload()
elapsed_ns = time_ns() - start_ns

report = Dict(
    "meta" => benchmark_metadata(run_tag="cold_start"),
    "policy" => "cold start / TTFx metric only; not for Python-vs-Julia runtime comparison",
    "first_call_ns" => elapsed_ns,
)

mkpath(joinpath(@__DIR__, "..", "..", "reports"))
out = joinpath(@__DIR__, "..", "..", "reports", "julia_cold_start.json")
open(out, "w") do io
    JSON3.write(io, report)
end
println(report)
