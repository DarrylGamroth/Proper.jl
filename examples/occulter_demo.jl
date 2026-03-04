using proper
using Plots
include(joinpath(@__DIR__, "_shared.jl"))

psf, sampling = run_simple_case()
println("occulter_demo: sampling = ", sampling)
plot_psf(psf; title="occulter_demo")
