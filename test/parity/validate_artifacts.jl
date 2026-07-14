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
    @test metadata["schema_version"] == 1

    baseline = metadata["baseline"]
    @test baseline["name"] == "PROPER Python"
    @test baseline["version"] == "3.3.4"
    @test startswith(String(baseline["source_url"]), "https://sourceforge.net/")
    @test occursin(r"^[0-9a-f]{64}$", String(baseline["source_snapshot_sha256"]))

    environment = metadata["environment"]
    for key in ("python", "python_executable", "numpy", "scipy", "astropy", "backend")
        @test haskey(environment, key)
        @test !isempty(String(environment[key]))
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
