using JSON3

include(joinpath(@__DIR__, "..", "..", "common", "metadata.jl"))
using .BenchMetadata

backend_name = lowercase(get(ENV, "PROPER_LATENCY_BACKEND", "unknown"))
backend = Symbol(backend_name)
precision_value = lowercase(get(ENV, "PROPER_LATENCY_PRECISION", "float64"))
precision_tag = if precision_value in ("float32", "fp32")
    "fp32"
elseif precision_value in ("float64", "fp64")
    "fp64"
else
    throw(ArgumentError("PROPER_LATENCY_PRECISION must be Float64/fp64 or Float32/fp32"))
end
reason = get(ENV, "PROPER_LATENCY_SKIP_REASON", "backend unavailable")
report_path = get(
    ENV,
    "PROPER_LATENCY_REPORT",
    joinpath(
        @__DIR__,
        "..",
        "..",
        "reports",
        "julia_prepared_latency_$(backend_name)_$(precision_tag).json",
    ),
)

report = Dict{String,Any}(
    "meta" => merge(
        benchmark_metadata(
            run_tag="prepared_latency_$(backend_name)_$(precision_tag)",
            backend=backend,
        ),
        Dict{String,Any}("available" => false),
    ),
    "policy" => "prepared-call latency benchmark skipped before measurement",
    "reason" => reason,
)

mkpath(dirname(report_path))
open(report_path, "w") do io
    JSON3.pretty(io, report)
end
println(JSON3.write(report))
