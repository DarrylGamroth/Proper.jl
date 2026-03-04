using proper
using Plots
include(joinpath(@__DIR__, "_shared.jl"))

psf, sampling = run_simple_case()
println("coronagraph: sampling = ", sampling)
plot_psf(psf; title="coronagraph")
