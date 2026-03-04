using proper
using Plots
include(joinpath(@__DIR__, "_shared.jl"))

psf, sampling = run_simple_case()
println("example_system: sampling = ", sampling)
plot_psf(psf; title="example_system")
