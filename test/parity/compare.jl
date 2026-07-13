using DelimitedFiles
using JSON3
using LinearAlgebra
using Proper
using TOML

include(joinpath(@__DIR__, "gates.jl"))

function run_simple_case()
    wf = prop_begin(2.4, 0.55e-6, 256; beam_diam_fraction=0.5)
    prop_circular_aperture(wf, 0.6)
    prop_lens(wf, 20.0)
    prop_propagate(wf, 20.0)
    return prop_end(wf)
end

base = joinpath(@__DIR__, "baseline", "python334")
psf_py = readdlm(joinpath(base, "simple_case_psf.csv"), ',')
meta = JSON3.read(read(joinpath(base, "simple_case_meta.json"), String))
psf_jl, sampling = run_simple_case()

expected_shape = Tuple(Int.(meta["shape"]))
for (label, actual_shape) in (("Julia output", size(psf_jl)), ("Python CSV", size(psf_py)))
    shape_failure = parity_shape_failure(actual_shape, expected_shape)
    isnothing(shape_failure) || error("Simple-case parity $label $shape_failure")
end

rel_l2 = norm(psf_jl .- psf_py) / max(norm(psf_py), eps())
expected_sampling = Float64(meta["sampling"])
sampling_relerr = abs(sampling - expected_sampling) / max(abs(expected_sampling), eps())

case_config = TOML.parsefile(joinpath(@__DIR__, "cases", "simple_case.toml"))
threshold_config = case_config["thresholds"]
metrics = Dict(
    "relative_l2" => Float64(rel_l2),
    "sampling_relerr" => Float64(sampling_relerr),
)
thresholds = Dict(
    "relative_l2" => Float64(threshold_config["relative_l2_max"]),
    "sampling_relerr" => Float64(threshold_config["sampling_relerr_max"]),
)
failures = parity_threshold_failures(metrics, thresholds)
report = Dict(
    "case" => "simple_case",
    "relative_l2" => rel_l2,
    "sampling" => sampling,
    "expected_sampling" => expected_sampling,
    "sampling_relerr" => sampling_relerr,
    "thresholds" => thresholds,
    "failures" => failures,
    "pass" => isempty(failures),
)

mkpath(joinpath(@__DIR__, "reports"))
open(joinpath(@__DIR__, "reports", "simple_case_report.json"), "w") do io
    JSON3.write(io, report)
end

println(report)
isempty(failures) || error("Simple-case parity threshold check failed: $(join(failures, "; "))")
