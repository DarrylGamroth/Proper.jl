using proper
using Random
using Plots

function run_simple_case(; diam=2.4, lambda_um=0.55, gridsize=256, beam_frac=0.5, central_obscuration=0.0, lens_fl=20.0, dz=20.0)
    wf = prop_begin(diam, lambda_um * 1e-6, gridsize; beam_diam_fraction=beam_frac)
    prop_circular_aperture(wf, diam * beam_frac / 2)
    if central_obscuration > 0
        prop_circular_obscuration(wf, central_obscuration * diam * beam_frac / 2)
    end
    prop_lens(wf, lens_fl)
    prop_propagate(wf, dz)
    return prop_end(wf)
end

function plot_psf(psf; title="PSF", logscale=true)
    img = logscale ? log10.(psf .+ eps(eltype(psf))) : psf
    plt = heatmap(img, aspect_ratio=:equal, colorbar=true, title=title)
    display(plt)
    return plt
end

function apply_seeded_psd!(wf; amp=1e-9, seed=42)
    rng = MersenneTwister(seed)
    prop_psd_errormap(wf, amp, 2.0, 3.0; rng=rng)
    return wf
end
