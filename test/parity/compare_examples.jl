using JSON3
using LinearAlgebra
using Random
using Proper

include(joinpath(@__DIR__, "..", "example_loader.jl"))
include(joinpath(@__DIR__, "gates.jl"))

const EXDIR = joinpath(@__DIR__, "..", "..", "examples")
const EXAMPLE_MODULES = Dict{Symbol,Module}(
    :simple_prescription => load_example_module(joinpath(EXDIR, "simple_prescription.jl")),
    :simple_telescope => load_example_module(joinpath(EXDIR, "simple_telescope.jl")),
    :hubble_simple => load_example_module(joinpath(EXDIR, "hubble_simple.jl")),
    :microscope => load_example_module(joinpath(EXDIR, "microscope.jl")),
    :example_system => load_example_module(joinpath(EXDIR, "example_system.jl")),
    :psdtest => load_example_module(joinpath(EXDIR, "psdtest.jl")),
    :talbot => load_example_module(joinpath(EXDIR, "talbot.jl")),
    :talbot_correct => load_example_module(joinpath(EXDIR, "talbot_correct.jl")),
    :run_occulter => load_example_module(joinpath(EXDIR, "run_occulter.jl")),
    :run_coronagraph => load_example_module(joinpath(EXDIR, "run_coronagraph.jl")),
    :run_coronagraph_dm => load_example_module(joinpath(EXDIR, "run_coronagraph_dm.jl")),
    :multi_example => load_example_module(joinpath(EXDIR, "multi_example.jl")),
    :testmulti1 => load_example_module(joinpath(EXDIR, "testmulti1.jl")),
    :testmulti2 => load_example_module(joinpath(EXDIR, "testmulti2.jl")),
)

function probe_indices(a::AbstractMatrix)
    ny, nx = size(a)
    return (
        (ny ÷ 5 + 1, nx ÷ 3 + 1),
        (ny ÷ 3 + 1, 2 * nx ÷ 3 + 1),
        (3 * ny ÷ 5 + 1, nx ÷ 4 + 1),
        (4 * ny ÷ 5 + 1, 3 * nx ÷ 4 + 1),
    )
end

function summarize(a, sampling)
    probes = [a[index...] for index in probe_indices(a)]
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
            "probes_re" => Float64[real(value) for value in probes],
            "probes_im" => Float64[imag(value) for value in probes],
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
            "probes" => Float64[probes...],
        )
    end
end

Random.seed!(12345)

testmulti1_psf = getfield(EXAMPLE_MODULES[:testmulti1], :testmulti1)()
testmulti2_fields, testmulti2_samplings = getfield(EXAMPLE_MODULES[:testmulti2], :testmulti2)()

cases = Dict{String,Tuple{Any,Any}}(
    "simple_prescription" => getfield(EXAMPLE_MODULES[:simple_prescription], :simple_prescription)(0.55e-6, 256),
    "simple_telescope" => getfield(EXAMPLE_MODULES[:simple_telescope], :simple_telescope)(0.55e-6, 256),
    "hubble_simple" => getfield(EXAMPLE_MODULES[:hubble_simple], :hubble_simple)(0.55e-6, 256, Dict("delta_sec" => 0.0)),
    "microscope" => getfield(EXAMPLE_MODULES[:microscope], :microscope)(0.55e-6, 256, Dict("focus_offset" => 0.0)),
    "example_system" => getfield(EXAMPLE_MODULES[:example_system], :example_system)(0.55e-6, 256),
    "psdtest" => getfield(EXAMPLE_MODULES[:psdtest], :psdtest)(0.55e-6, 256, Dict("usepsdmap" => true)),
    "talbot" => getfield(EXAMPLE_MODULES[:talbot], :talbot)(0.5e-6, 128, Dict("diam" => 0.1, "period" => 0.04, "dist" => 0.0)),
    "talbot_correct" => getfield(EXAMPLE_MODULES[:talbot_correct], :talbot_correct)(0.5e-6, 128, Dict("diam" => 0.1, "period" => 0.04, "dist" => 0.0)),
    "run_occulter" => getfield(EXAMPLE_MODULES[:run_occulter], :run_occulter)(0.55e-6, 256, Dict("occulter_type" => "GAUSSIAN")),
    "run_coronagraph" => getfield(EXAMPLE_MODULES[:run_coronagraph], :run_coronagraph)(0.55e-6, 256, Dict("use_errors" => false, "occulter_type" => "GAUSSIAN")),
    "run_coronagraph_dm" => getfield(EXAMPLE_MODULES[:run_coronagraph_dm], :run_coronagraph_dm)(0.55e-6, 256, Dict("use_errors" => false, "use_dm" => false, "occulter_type" => "GAUSSIAN")),
    "multi_example" => getfield(EXAMPLE_MODULES[:multi_example], :multi_example)(0.55e-6, 256, Dict("use_dm" => false, "dm" => zeros(48, 48))),
    "testmulti1" => (testmulti1_psf, 1.5e-6),
    "testmulti2_pattern_1" => (testmulti2_fields[:, :, 1], testmulti2_samplings[1]),
    "testmulti2_pattern_2" => (testmulti2_fields[:, :, 2], testmulti2_samplings[2]),
    "testmulti2_pattern_3" => (testmulti2_fields[:, :, 3], testmulti2_samplings[3]),
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


function vector_max_absdiff(actual, expected)
    length(actual) == length(expected) || return Inf
    return maximum(abs(Float64(actual[i]) - Float64(expected[i])) for i in eachindex(actual))
end


function vector_max_relerr_floor(actual, expected; floor::Float64=1e-12)
    length(actual) == length(expected) || return Inf
    return maximum(
        relerr_floor(Float64(actual[i]), Float64(expected[i]); floor)
        for i in eachindex(actual)
    )
end

report = Dict{String,Any}()
structural_failures = Dict{String,Vector{String}}()
for (name, j) in jl
    p = py[name]
    shape_failure = parity_shape_failure(Tuple(Int.(j["shape"])), Tuple(Int.(p["shape"])))
    if !isnothing(shape_failure)
        structural_failures[name] = [shape_failure]
    end
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
        m["center_absdiff"] = abs(j["center"] - Float64(p["center"]))
        m["center_relerr_floor"] = relerr_floor(j["center"], Float64(p["center"]))
        m["probe_absdiff"] = vector_max_absdiff(j["probes"], p["probes"])
        m["probe_relerr_floor"] = vector_max_relerr_floor(j["probes"], p["probes"])
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
        m["center_re_absdiff"] = abs(j["center_re"] - Float64(p["center_re"]))
        m["center_im_absdiff"] = abs(j["center_im"] - Float64(p["center_im"]))
        m["center_re_relerr_floor"] = relerr_floor(j["center_re"], Float64(p["center_re"]))
        m["center_im_relerr_floor"] = relerr_floor(j["center_im"], Float64(p["center_im"]))
        m["probe_re_absdiff"] = vector_max_absdiff(j["probes_re"], p["probes_re"])
        m["probe_im_absdiff"] = vector_max_absdiff(j["probes_im"], p["probes_im"])
        m["probe_re_relerr_floor"] = vector_max_relerr_floor(j["probes_re"], p["probes_re"])
        m["probe_im_relerr_floor"] = vector_max_relerr_floor(j["probes_im"], p["probes_im"])
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
    metric_thresholds = Dict{String,Float64}()
    for (kmax, vmax) in th
        metric = endswith(kmax, "_max") ? kmax[1:(end - 4)] : kmax
        metric_thresholds[metric] = vmax
    end
    bad = copy(get(structural_failures, name, String[]))
    append!(bad, parity_threshold_failures(m, metric_thresholds))
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
