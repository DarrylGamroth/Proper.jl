using proper
using Plots
include(joinpath(@__DIR__, "_shared.jl"))

psf, sampling = run_simple_case()
println("talbot_demo: sampling = ", sampling)
plot_psf(psf; title="talbot_demo")
