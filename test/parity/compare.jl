using DelimitedFiles
using JSON3
using LinearAlgebra
using Proper
using Statistics
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

@inline json_complex_scalar(payload) = complex(
    Float64(payload["real"]),
    Float64(payload["imag"]),
)

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

function compare_carrier_phase(base)
    baseline = JSON3.read(read(joinpath(base, "carrier_phase.json"), String))
    config = TOML.parsefile(joinpath(@__DIR__, "cases", "carrier_phase.toml"))
    field_threshold = Float64(config["thresholds"]["field_max_abs"])
    intensity_threshold = Float64(config["thresholds"]["intensity_max_abs"])
    n = Int(baseline["gridsize"])
    λ = Float64(baseline["wavelength_m"])

    function propagated(distance, carrier_phase)
        wf = prop_begin(1.0, λ, n)
        ctx = RunContext(
            typeof(wf.field),
            wf.workspace;
            carrier_phase=carrier_phase,
        )
        prop_ptp(wf, distance, ctx)
        return wf.field
    end

    quarter_disabled = propagated(λ / 4, EnvelopeOnly())
    quarter_enabled = propagated(λ / 4, TrackCarrierPhase())
    arm1 = propagated(λ, TrackCarrierPhase())
    arm2 = propagated(3λ / 2, TrackCarrierPhase())
    arm_sum = arm1 .+ arm2

    actual = Dict(
        "quarter_disabled_mean" => mean(quarter_disabled),
        "quarter_enabled_mean" => mean(quarter_enabled),
        "arm_sum_mean" => mean(arm_sum),
        "arm_mean_intensity" => mean(abs2, arm_sum),
    )
    outputs = baseline["outputs"]
    expected = Dict(
        "quarter_disabled_mean" => json_complex_scalar(outputs["quarter_disabled_mean"]),
        "quarter_enabled_mean" => json_complex_scalar(outputs["quarter_enabled_mean"]),
        "arm_sum_mean" => json_complex_scalar(outputs["arm_sum_mean"]),
        "arm_mean_intensity" => Float64(outputs["arm_mean_intensity"]),
    )

    metric_type = @NamedTuple{max_abs::Float64, threshold::Float64}
    metrics = Dict{String,metric_type}()
    failures = String[]
    for label in keys(actual)
        threshold = label == "arm_mean_intensity" ? intensity_threshold : field_threshold
        max_abs = Float64(abs(actual[label] - expected[label]))
        metrics[label] = (max_abs=max_abs, threshold=threshold)
        max_abs <= threshold || push!(failures, "$label max_abs=$max_abs exceeds $threshold")
    end

    report = (
        case="carrier_phase",
        metrics=metrics,
        failures=failures,
        pass=isempty(failures),
    )
    mkpath(joinpath(@__DIR__, "reports"))
    open(joinpath(@__DIR__, "reports", "carrier_phase_report.json"), "w") do io
        JSON3.write(io, report)
    end
    println(report)
    isempty(failures) || error("Carrier-phase parity failed: $(join(failures, "; "))")
    return report
end

function compare_dm_tilt_projection(base)
    baseline = JSON3.read(read(joinpath(base, "dm_tilt_projection.json"), String))
    thresholds = TOML.parsefile(joinpath(@__DIR__, "cases", "dm_tilt_projection.toml"))["thresholds"]
    config = baseline["config"]
    dm_surface = json_real_matrix(baseline["inputs"]["dm_surface_m"])
    n = Int(config["grid_size"])
    expected_shape = Tuple(Int.(baseline["shape"]))
    expected_shape == (n, n) || error("DM tilt baseline shape/config mismatch")

    function wavefront()
        return prop_begin(
            Float64(config["beam_diameter_m"]),
            Float64(config["wavelength_m"]),
            n;
            beam_diam_fraction=Float64(config["beam_diam_fraction"]),
        )
    end

    kwargs = (
        XTILT=Float64(config["xtilt_deg"]),
        YTILT=Float64(config["ytilt_deg"]),
        ZTILT=Float64(config["ztilt_deg"]),
    )
    dm_xc = Float64(config["dm_xc"])
    dm_yc = Float64(config["dm_yc"])
    spacing = Float64(config["spacing_m"])
    wf_map = wavefront()
    actual_dmap = prop_dm(
        wf_map,
        dm_surface,
        dm_xc,
        dm_yc,
        spacing;
        NO_APPLY=true,
        kwargs...,
    )
    wf_zero_tilt = wavefront()
    actual_dmap_zero_tilt = prop_dm(
        wf_zero_tilt,
        dm_surface,
        dm_xc,
        dm_yc,
        spacing;
        NO_APPLY=true,
    )
    wf_apply = wavefront()
    prop_dm(wf_apply, dm_surface, dm_xc, dm_yc, spacing; kwargs...)
    actual_wavefront = prop_get_wavefront(wf_apply)

    outputs = baseline["outputs"]
    expected_dmap_zero_tilt = json_real_matrix(outputs["dmap_zero_tilt_m"])
    expected_dmap = json_real_matrix(outputs["dmap_m"])
    expected_wavefront = json_complex_matrix(outputs["centered_wavefront"])
    metric_type = @NamedTuple{
        max_abs::Float64,
        worst_index::NTuple{2,Int},
        threshold::Float64,
    }
    metrics = Dict{String,metric_type}()
    failures = String[]
    record_pixelwise_comparison!(
        metrics,
        failures,
        "dmap_zero_tilt_m",
        actual_dmap_zero_tilt,
        expected_dmap_zero_tilt,
        Float64(thresholds["dmap_max_abs"]),
    )
    record_pixelwise_comparison!(
        metrics,
        failures,
        "dmap_m",
        actual_dmap,
        expected_dmap,
        Float64(thresholds["dmap_max_abs"]),
    )
    record_pixelwise_comparison!(
        metrics,
        failures,
        "centered_wavefront",
        actual_wavefront,
        expected_wavefront,
        Float64(thresholds["field_max_abs"]),
    )

    report = (
        case="dm_tilt_projection",
        shape=expected_shape,
        metrics=metrics,
        failures=failures,
        pass=isempty(failures),
    )
    mkpath(joinpath(@__DIR__, "reports"))
    open(joinpath(@__DIR__, "reports", "dm_tilt_projection_report.json"), "w") do io
        JSON3.write(io, report)
    end
    println(report)
    isempty(failures) || error("Tilted-DM parity failed: $(join(failures, "; "))")
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
compare_carrier_phase(base)
compare_dm_tilt_projection(base)
