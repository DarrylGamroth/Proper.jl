using proper
using Plots
include(joinpath(@__DIR__, "_shared.jl"))

wf = prop_begin(2.4, 0.55e-6, 256; beam_diam_fraction=0.5)
prop_circular_aperture(wf, 0.6)
apply_seeded_psd!(wf; amp=5e-10, seed=123)
psf, sampling = prop_end(wf)
println("psdtest: sampling = ", sampling)
plot_psf(psf; title="psdtest")
