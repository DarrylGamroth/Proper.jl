using proper
using Plots
include(joinpath(@__DIR__, "_shared.jl"))

psf, sampling = run_simple_case()
println("hubble_simple: sampling = ", sampling)
plot_psf(psf; title="hubble_simple")
