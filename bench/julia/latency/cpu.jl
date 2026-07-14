include(joinpath(@__DIR__, "..", "..", "common", "prepared_latency_report.jl"))

run_prepared_latency_report(:cpu, cpu_latency_field, () -> nothing)
