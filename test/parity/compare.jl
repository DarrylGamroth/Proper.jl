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

function json_real_matrix(rows)
    ny = length(rows)
    ny > 0 || error("Parity matrix cannot be empty")
    nx = length(rows[1])
    all(length(row) == nx for row in rows) || error("Parity matrix rows must have equal lengths")
    return [Float64(rows[i][j]) for i in 1:ny, j in 1:nx]
end

function json_complex_matrix(payload)
    real_part = json_real_matrix(payload["real"])
    imag_part = json_real_matrix(payload["imag"])
    size(real_part) == size(imag_part) || error("Complex parity components must have equal shapes")
    return complex.(real_part, imag_part)
end

function record_pixelwise_comparison!(metrics, failures, label, actual, expected, max_abs)
    if size(actual) != size(expected)
        push!(failures, "$label shape mismatch: got $(size(actual)), expected $(size(expected))")
        return
    end

    errors = abs.(actual .- expected)
    if !all(isfinite, errors)
        push!(failures, "$label contains a non-finite pixel error")
        return
    end

    worst_error, worst_index = findmax(errors)
    metrics[label] = (
        max_abs=Float64(worst_error),
        worst_index=Tuple(worst_index),
        threshold=Float64(max_abs),
    )
    worst_error <= max_abs || push!(
        failures,
        "$label max_abs=$(Float64(worst_error)) exceeds $(Float64(max_abs)) at $(Tuple(worst_index))",
    )
    return
end

function compare_wavefront_accessors_even(base)
    baseline_path = joinpath(base, "wavefront_accessors_even.json")
    baseline = JSON3.read(read(baseline_path, String))
    config = TOML.parsefile(joinpath(@__DIR__, "cases", "wavefront_accessors_even.toml"))
    thresholds = config["thresholds"]
    field_max_abs = Float64(thresholds["field_max_abs"])
    derived_max_abs = Float64(thresholds["derived_max_abs"])

    expected_shape = Tuple(Int.(baseline["shape"]))
    expected_shape == (Int(config["gridsize"]), Int(config["gridsize"])) ||
        error("wavefront accessor case/config shape mismatch")

    inputs = baseline["inputs"]
    raw = json_complex_matrix(inputs["raw_field"])
    centered_addition = json_complex_matrix(inputs["centered_addition"])
    scalar_payload = inputs["scalar_addition"]
    scalar_addition = complex(
        Float64(scalar_payload["real"]),
        Float64(scalar_payload["imag"]),
    )
    for (label, matrix) in (("raw field", raw), ("centered addition", centered_addition))
        size(matrix) == expected_shape || error("$label shape mismatch in Python baseline")
    end

    wf_accessors = Proper.WaveFront(copy(raw), 500e-9, 1e-3, 0.0, 1.0)
    actual_wavefront = prop_get_wavefront(wf_accessors)
    actual_amplitude = prop_get_amplitude(wf_accessors)
    actual_phase = prop_get_phase(wf_accessors)

    wf_scalar = Proper.WaveFront(copy(raw), 500e-9, 1e-3, 0.0, 1.0)
    prop_add_wavefront(wf_scalar, scalar_addition)
    actual_scalar_add = copy(wf_scalar.field)

    wf_matrix = Proper.WaveFront(copy(raw), 500e-9, 1e-3, 0.0, 1.0)
    prop_add_wavefront(wf_matrix, centered_addition)
    actual_matrix_add = copy(wf_matrix.field)

    outputs = baseline["outputs"]
    expected_wavefront = json_complex_matrix(outputs["get_wavefront"])
    expected_amplitude = json_real_matrix(outputs["get_amplitude"])
    expected_phase = json_real_matrix(outputs["get_phase"])
    expected_scalar_add = json_complex_matrix(outputs["after_scalar_add"])
    expected_matrix_add = json_complex_matrix(outputs["after_matrix_add"])

    metric_type = @NamedTuple{
        max_abs::Float64,
        worst_index::NTuple{2,Int},
        threshold::Float64,
    }
    metrics = Dict{String,metric_type}()
    failures = String[]
    for (label, actual, expected, threshold) in (
        ("get_wavefront", actual_wavefront, expected_wavefront, field_max_abs),
        ("get_amplitude", actual_amplitude, expected_amplitude, derived_max_abs),
        ("get_phase", actual_phase, expected_phase, derived_max_abs),
        ("after_scalar_add", actual_scalar_add, expected_scalar_add, field_max_abs),
        ("after_matrix_add", actual_matrix_add, expected_matrix_add, field_max_abs),
    )
        record_pixelwise_comparison!(metrics, failures, label, actual, expected, threshold)
    end

    report = (
        case="wavefront_accessors_even",
        shape=expected_shape,
        metrics=metrics,
        failures=failures,
        pass=isempty(failures),
    )
    mkpath(joinpath(@__DIR__, "reports"))
    open(joinpath(@__DIR__, "reports", "wavefront_accessors_even_report.json"), "w") do io
        JSON3.write(io, report)
    end
    println(report)
    isempty(failures) || error(
        "Wavefront accessor/add parity failed: $(join(failures, "; "))",
    )
    return report
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

compare_wavefront_accessors_even(base)
