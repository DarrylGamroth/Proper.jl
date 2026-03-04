using proper
using Plots
include(joinpath(@__DIR__, "_shared.jl"))

psf, sampling = run_simple_case()
println("simple_prescription: sampling = ", sampling)
plot_psf(psf; title="simple_prescription")
