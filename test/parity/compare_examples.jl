using JSON3
using LinearAlgebra
using Random
using proper

include(joinpath(@__DIR__, "..", "..", "examples", "simple_prescription.jl"))
include(joinpath(@__DIR__, "..", "..", "examples", "simple_telescope.jl"))
include(joinpath(@__DIR__, "..", "..", "examples", "hubble_simple.jl"))
include(joinpath(@__DIR__, "..", "..", "examples", "microscope.jl"))
include(joinpath(@__DIR__, "..", "..", "examples", "example_system.jl"))
include(joinpath(@__DIR__, "..", "..", "examples", "psdtest.jl"))
include(joinpath(@__DIR__, "..", "..", "examples", "talbot.jl"))
include(joinpath(@__DIR__, "..", "..", "examples", "talbot_correct.jl"))
include(joinpath(@__DIR__, "..", "..", "examples", "run_occulter.jl"))
include(joinpath(@__DIR__, "..", "..", "examples", "run_coronagraph.jl"))
include(joinpath(@__DIR__, "..", "..", "examples", "run_coronagraph_dm.jl"))
include(joinpath(@__DIR__, "..", "..", "examples", "multi_example.jl"))

function summarize(a, sampling)
    if eltype(a) <: Complex
        c = a[size(a, 1) ÷ 2 + 1, size(a, 2) ÷ 2 + 1]
        return Dict(
            "shape" => [size(a, 1), size(a, 2)],
            "sampling" => float(sampling),
            "sum_re" => float(sum(real, a)),
            "sum_im" => float(sum(imag, a)),
            "norm" => float(norm(a)),
            "center_re" => float(real(c)),
            "center_im" => float(imag(c)),
        )
    else
        c = a[size(a, 1) ÷ 2 + 1, size(a, 2) ÷ 2 + 1]
        return Dict(
            "shape" => [size(a, 1), size(a, 2)],
            "sampling" => float(sampling),
            "sum" => float(sum(a)),
            "norm" => float(norm(a)),
            "max" => float(maximum(a)),
            "center" => float(c),
        )
    end
end

Random.seed!(12345)

cases = Dict{String,Tuple{Any,Any}}(
    "simple_prescription" => simple_prescription(0.55e-6, 256),
    "simple_telescope" => simple_telescope(0.55e-6, 256),
    "hubble_simple" => hubble_simple(0.55e-6, 256, Dict("delta_sec" => 0.0)),
    "microscope" => microscope(0.55e-6, 256, Dict("focus_offset" => 0.0)),
    "example_system" => example_system(0.55e-6, 256),
    "psdtest" => psdtest(0.55e-6, 256, Dict("usepsdmap" => true)),
    "talbot" => talbot(0.5e-6, 128, Dict("diam" => 0.1, "period" => 0.04, "dist" => 0.0)),
    "talbot_correct" => talbot_correct(0.5e-6, 128, Dict("diam" => 0.1, "period" => 0.04, "dist" => 0.0)),
    "run_occulter" => run_occulter(0.55e-6, 256, Dict("occulter_type" => "GAUSSIAN")),
    "run_coronagraph" => run_coronagraph(0.55e-6, 256, Dict("use_errors" => false, "occulter_type" => "GAUSSIAN")),
    "run_coronagraph_dm" => run_coronagraph_dm(0.55e-6, 256, Dict("use_errors" => false, "use_dm" => false, "occulter_type" => "GAUSSIAN")),
    "multi_example" => multi_example(0.55e-6, 256, Dict("use_dm" => false, "dm" => zeros(48, 48))),
)

jl = Dict(name => summarize(psf, samp) for (name, (psf, samp)) in cases)

base = joinpath(@__DIR__, "baseline", "python334")
py = JSON3.read(read(joinpath(base, "example_metrics.json"), String))

function relerr(a::Real, b::Real)
    den = max(abs(b), eps())
    return abs(a - b) / den
end

function relerr_floor(a::Real, b::Real; floor::Float64=1e-12)
    den = max(abs(b), floor)
    return abs(a - b) / den
end

report = Dict{String,Any}()
for (name, j) in jl
    p = py[name]
    m = Dict{String,Any}()
    m["sampling_relerr"] = relerr(j["sampling"], Float64(p["sampling"]))
    if haskey(j, "sum")
        m["sum_relerr"] = relerr(j["sum"], Float64(p["sum"]))
        m["max_relerr"] = relerr(j["max"], Float64(p["max"]))
        m["norm_relerr"] = relerr(j["norm"], Float64(p["norm"]))
        m["sum_absdiff"] = abs(j["sum"] - Float64(p["sum"]))
        m["max_absdiff"] = abs(j["max"] - Float64(p["max"]))
        m["norm_absdiff"] = abs(j["norm"] - Float64(p["norm"]))
        m["sum_relerr_floor"] = relerr_floor(j["sum"], Float64(p["sum"]))
        m["max_relerr_floor"] = relerr_floor(j["max"], Float64(p["max"]))
        m["norm_relerr_floor"] = relerr_floor(j["norm"], Float64(p["norm"]))
    else
        m["sum_re_relerr"] = relerr(j["sum_re"], Float64(p["sum_re"]))
        m["sum_im_relerr"] = relerr(j["sum_im"], Float64(p["sum_im"]))
        m["norm_relerr"] = relerr(j["norm"], Float64(p["norm"]))
        m["sum_re_absdiff"] = abs(j["sum_re"] - Float64(p["sum_re"]))
        m["sum_im_absdiff"] = abs(j["sum_im"] - Float64(p["sum_im"]))
        m["norm_absdiff"] = abs(j["norm"] - Float64(p["norm"]))
        m["sum_re_relerr_floor"] = relerr_floor(j["sum_re"], Float64(p["sum_re"]))
        m["sum_im_relerr_floor"] = relerr_floor(j["sum_im"], Float64(p["sum_im"]))
        m["norm_relerr_floor"] = relerr_floor(j["norm"], Float64(p["norm"]))
    end
    report[name] = m
end

mkpath(joinpath(@__DIR__, "reports"))
open(joinpath(@__DIR__, "reports", "example_metrics_report.json"), "w") do io
    JSON3.write(io, report)
end

thresholds = JSON3.read(read(joinpath(@__DIR__, "thresholds", "example_metrics_thresholds.json"), String))
default_real = Dict{String,Float64}(String(k) => Float64(v) for (k, v) in pairs(thresholds["default_real"]))
default_complex = Dict{String,Float64}(String(k) => Float64(v) for (k, v) in pairs(thresholds["default_complex"]))
overrides = Dict{String,Dict{String,Float64}}()
for (k, v) in pairs(thresholds["overrides"])
    overrides[String(k)] = Dict{String,Float64}(String(kk) => Float64(vv) for (kk, vv) in pairs(v))
end

function merged_thresholds(name::String, is_complex::Bool)
    base = is_complex ? copy(default_complex) : copy(default_real)
    if haskey(overrides, name)
        merge!(base, overrides[name])
    end
    return base
end

failures = Dict{String,Any}()
for (name, m_any) in report
    m = Dict{String,Float64}(String(k) => Float64(v) for (k, v) in m_any)
    is_complex = haskey(m, "sum_re_relerr")
    th = merged_thresholds(name, is_complex)
    bad = String[]
    for (kmax, vmax) in th
        metric = endswith(kmax, "_max") ? kmax[1:(end - 4)] : kmax
        if !haskey(m, metric)
            push!(bad, "missing metric $metric")
            continue
        end
        if m[metric] > vmax
            push!(bad, "$metric=$(m[metric]) > $vmax")
        end
    end
    if !isempty(bad)
        failures[name] = bad
    end
end

summary = Dict(
    "pass" => isempty(failures),
    "num_cases" => length(report),
    "num_failed" => length(failures),
    "failures" => failures,
    "report_file" => joinpath(@__DIR__, "reports", "example_metrics_report.json"),
)
open(joinpath(@__DIR__, "reports", "example_metrics_threshold_summary.json"), "w") do io
    JSON3.write(io, summary)
end

println(report)
println(summary)
isempty(failures) || error("Example parity threshold check failed for $(length(failures)) case(s).")
