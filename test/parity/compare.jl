using DelimitedFiles
using JSON3
using LinearAlgebra
using Proper

function run_simple_case()
    wf = prop_begin(2.4, 0.55e-6, 256; beam_diam_fraction=0.5)
    prop_circular_aperture(wf, 0.6)
    prop_lens(wf, 20.0)
    prop_propagate(wf, 20.0)
    return prop_end(wf)
end

base = joinpath(@__DIR__, "baseline", "python334")
psf_py = readdlm(joinpath(base, "simple_case_psf.csv"), ',')
psf_jl, sampling = run_simple_case()

rel_l2 = norm(psf_jl .- psf_py) / max(norm(psf_py), eps())
report = Dict(
    "case" => "simple_case",
    "relative_l2" => rel_l2,
    "sampling" => sampling,
)

mkpath(joinpath(@__DIR__, "reports"))
open(joinpath(@__DIR__, "reports", "simple_case_report.json"), "w") do io
    JSON3.write(io, report)
end

println(report)
