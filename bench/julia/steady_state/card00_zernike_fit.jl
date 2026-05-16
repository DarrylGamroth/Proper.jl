using Proper

include(joinpath(@__DIR__, "..", "..", "common", "metadata.jl"))
include(joinpath(@__DIR__, "..", "..", "cards", "card00", "benchmark_helpers.jl"))
include(joinpath(@__DIR__, "..", "..", "cards", "card00", "reporting.jl"))

using .BenchMetadata
using .Card00BenchmarkHelpers
using .Card00Reporting

function bench_zernike_fit_case(grid_n::Integer, nzer::Integer, samples::Integer)
    wavefront = zernike_fit_wavefront(grid_n, max(Int(nzer), maximum(CARD00_ZERNIKE_NTERMS)))
    pupil = circular_fit_pupil(grid_n)

    coeff_workload() = begin
        prop_fit_zernikes(wavefront, pupil.mask, pupil.radius, nzer; xc=pupil.xc, yc=pupil.yc)
        nothing
    end
    fitmap_workload() = begin
        prop_fit_zernikes(wavefront, pupil.mask, pupil.radius, nzer; xc=pupil.xc, yc=pupil.yc, fit=true)
        nothing
    end

    return Dict(
        "grid" => Int(grid_n),
        "nzer" => Int(nzer),
        "sample_count" => count(!iszero, pupil.mask),
        "coefficients" => trial_stats(measure_samples(coeff_workload, samples)),
        "coefficients_and_fit_map" => trial_stats(measure_samples(fitmap_workload, samples)),
    )
end

function main()
    samples = parse(Int, String(arg_value("--samples", "5")))
    grids = parse_int_list(arg_value("--grids", nothing), CARD00_ZERNIKE_FIT_GRIDS)
    nzer_list = parse_int_list(arg_value("--nzer", nothing), CARD00_ZERNIKE_NTERMS)

    cases = [
        bench_zernike_fit_case(grid_n, nzer, samples)
        for grid_n in grids for nzer in nzer_list
    ]
    report = Dict(
        "meta" => merge(
            benchmark_metadata(run_tag="card_00_zernike_fit", backend=:cpu),
            Dict("card" => "00", "benchmark" => "zernike_fit"),
        ),
        "policy" => "Card 00 steady-state CPU prop_fit_zernikes benchmarks via BenchmarkTools with evals=1 after warmup; TTFx excluded. Wavefront map, pupil samples, and radius are deterministic per grid.",
        "cases" => cases,
    )
    out = joinpath(@__DIR__, "..", "..", "reports", "card_00_zernike_fit.json")
    write_json_report(out, report)
end

main()
