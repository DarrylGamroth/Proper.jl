using proper
using Plots
using Statistics

function plot_psf(psf; title="PSF", power=0.25)
    img = psf .^ power
    plt = heatmap(img; aspect_ratio=:equal, colorbar=true, title=title)
    display(plt)
    return plt
end

function center_crop(a::AbstractMatrix, n::Integer)
    ny, nx = size(a)
    cy = ny ÷ 2 + 1
    cx = nx ÷ 2 + 1
    ry = (cy - n ÷ 2):(cy + (n - 1) ÷ 2)
    rx = (cx - n ÷ 2):(cx + (n - 1) ÷ 2)
    return @view a[ry, rx]
end

function maybe_save_plot(plt, path::AbstractString="")
    if !isempty(path)
        savefig(plt, path)
    end
    return nothing
end
