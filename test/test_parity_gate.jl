using Test

include(joinpath(@__DIR__, "parity", "gates.jl"))

@testset "Parity threshold gate" begin
    thresholds = Dict("relative_l2" => 1e-12, "sampling_relerr" => 1e-12)

    @test isempty(parity_threshold_failures(
        Dict("relative_l2" => 1e-15, "sampling_relerr" => 0.0),
        thresholds,
    ))
    @test isempty(parity_threshold_failures(
        Dict("relative_l2" => 1e-12, "sampling_relerr" => 1e-12),
        thresholds,
    ))

    exceeded = parity_threshold_failures(
        Dict("relative_l2" => 1e-3, "sampling_relerr" => 0.0),
        thresholds,
    )
    @test length(exceeded) == 1
    @test occursin("relative_l2", only(exceeded))

    missing = parity_threshold_failures(Dict("relative_l2" => 0.0), thresholds)
    @test length(missing) == 1
    @test occursin("missing metric sampling_relerr", only(missing))

    nonfinite = parity_threshold_failures(
        Dict("relative_l2" => NaN, "sampling_relerr" => Inf),
        thresholds,
    )
    @test length(nonfinite) == 2
    @test all(failure -> occursin("not finite", failure), nonfinite)

    invalid_thresholds = parity_threshold_failures(
        Dict("relative_l2" => 0.0, "sampling_relerr" => 0.0),
        Dict("relative_l2" => -1.0, "sampling_relerr" => Inf),
    )
    @test length(invalid_thresholds) == 2
    @test all(failure -> occursin("invalid threshold", failure), invalid_thresholds)

    @test isnothing(parity_shape_failure((256, 256), (256, 256)))
    @test occursin("shape mismatch", parity_shape_failure((1, 256), (256, 256)))
end
