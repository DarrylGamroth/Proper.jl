using JSON3
using SHA
using Test

const BASELINE_DIR = joinpath(@__DIR__, "baseline", "python334")
const REQUIRED_TOP_LEVEL_KEYS = (
    "schema_version",
    "baseline",
    "generator",
    "generated_at_utc",
    "environment",
    "native_kernels",
    "numeric_precision",
    "rng_seed",
    "cases",
    "artifacts",
)

function artifact_sha256(path::AbstractString)
    return bytes2hex(sha256(read(path)))
end

function validate_metadata(path::AbstractString)
    metadata = JSON3.read(read(path, String))
    @test all(key -> haskey(metadata, key), REQUIRED_TOP_LEVEL_KEYS)
    @test metadata["schema_version"] == 2

    baseline = metadata["baseline"]
    @test baseline["name"] == "PROPER Python"
    @test baseline["version"] == "3.3.4"
    @test startswith(String(baseline["source_url"]), "https://sourceforge.net/")
    expected_snapshot = "2f6f351715a49524f01aded1baedf7ef9c41bb40a3f738a4e7f481ed1daa7382"
    @test String(baseline["source_snapshot_sha256"]) == expected_snapshot

    bootstrap = baseline["reproducible_bootstrap"]
    @test bootstrap["seed_version"] == "3.3.5"
    @test startswith(String(bootstrap["seed_url"]), "https://sourceforge.net/")
    @test bootstrap["seed_archive_sha256"] == "5c25bc4ca80efb088990f1d6be231fe5583a806ff806de17ddc26026f2b23d87"
    @test bootstrap["expected_source_snapshot_sha256"] == expected_snapshot
    reconstruction_patch = normpath(joinpath(@__DIR__, "..", "..", String(bootstrap["reconstruction_patch"])))
    @test isfile(reconstruction_patch)
    @test artifact_sha256(reconstruction_patch) == String(bootstrap["reconstruction_patch_sha256"])

    environment = metadata["environment"]
    for key in ("python", "python_executable", "numpy", "scipy", "astropy", "backend")
        @test haskey(environment, key)
        @test !isempty(String(environment[key]))
    end
    @test environment["backend"] == "NumPy + upstream PROPER native C kernels"

    native_runtime_metadata = metadata["native_kernels"]
    @test native_runtime_metadata["required"] == true
    @test native_runtime_metadata["implementation"] == "upstream PROPER 3.3.4 C sources"
    @test !isempty(String(native_runtime_metadata["platform"]))
    @test !isempty(String(native_runtime_metadata["machine"]))
    @test occursin("threaded kernel assumes X", String(native_runtime_metadata["coordinate_grid_policy"]))
    compiler = native_runtime_metadata["compiler"]
    @test !isempty(compiler["command"])
    @test !isempty(String(compiler["version"]))
    kernels = native_runtime_metadata["kernels"]
    @test Set(String.(keys(kernels))) == Set(("cubic_conv", "cubic_conv_threaded", "szoom"))
    source_root = String(baseline["source_root"])
    for (name, kernel) in pairs(kernels)
        source = String(kernel["source"])
        @test startswith(source, "proper/")
        @test artifact_sha256(joinpath(source_root, source)) == String(kernel["source_sha256"])
        flags = Set(String.(kernel["compile_flags"]))
        @test "-ffp-contract=off" in flags
        @test "-fno-fast-math" in flags
        @test "-lm" in flags
        @test !isempty(String(kernel["required_symbol"]))
        if String(name) == "cubic_conv_threaded"
            @test "-pthread" in flags
        end
    end

    @test !isempty(String(metadata["generator"]))
    @test endswith(String(metadata["generated_at_utc"]), "Z")
    @test Set(String.(metadata["numeric_precision"])) == Set(("float64", "complex128"))

    cases = metadata["cases"]
    @test !isempty(cases)
    for (_, case) in pairs(cases)
        @test haskey(case, "wavelength_m")
        @test haskey(case, "grid_size")
        @test haskey(case, "config")
        @test Int(case["grid_size"]) > 0
    end

    artifacts = metadata["artifacts"]
    @test !isempty(artifacts)
    for (name, expected_hash) in pairs(artifacts)
        artifact_path = joinpath(dirname(path), String(name))
        @test isfile(artifact_path)
        @test artifact_sha256(artifact_path) == String(expected_hash)
    end
    return nothing
end

@testset "Python parity artifact provenance" begin
    validate_metadata(joinpath(BASELINE_DIR, "core_baseline_metadata.json"))
    validate_metadata(joinpath(BASELINE_DIR, "example_metrics_metadata.json"))
end
