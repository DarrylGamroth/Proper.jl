using Proper

function simple_prescription(wavelength::Real, gridsize::Integer)
    diam = 1.0
    focal_ratio = 15.0
    focal_length = diam * focal_ratio
    beam_ratio = 0.5

    wfo = prop_begin(diam, wavelength, gridsize; beam_diam_fraction=beam_ratio)
    prop_circular_aperture(wfo, diam / 2)
    prop_define_entrance(wfo)
    prop_lens(wfo, focal_length)
    prop_propagate(wfo, focal_length)
    return prop_end(wfo)
end

if abspath(PROGRAM_FILE) == @__FILE__
    using Plots

    psf, sampling = simple_prescription(0.55e-6, 256)
    println("simple_prescription: sampling = ", sampling)
    display(heatmap(psf .^ 0.25; aspect_ratio=:equal, color=:grays, title="simple_prescription"))
end
