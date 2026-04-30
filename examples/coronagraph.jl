using Proper
using Plots

abstract type AbstractOcculter end
struct GaussianOcculter <: AbstractOcculter end
struct SolidOcculter <: AbstractOcculter end
struct EighthOrderOcculter <: AbstractOcculter end

normalize_occulter(occulter::AbstractOcculter) = occulter

function normalize_occulter(occulter::Symbol)
    occulter in (:gaussian, :GaussianOcculter) && return GaussianOcculter()
    occulter in (:solid, :SolidOcculter) && return SolidOcculter()
    occulter in (:eighth_order, :EighthOrderOcculter) && return EighthOrderOcculter()
    throw(ArgumentError("Unknown occulter selector: $(occulter)"))
end

function normalize_occulter(occulter::AbstractString)
    normalized = lowercase(replace(occulter, '-' => '_', ' ' => '_'))
    normalized == "gaussian" && return GaussianOcculter()
    normalized == "solid" && return SolidOcculter()
    normalized in ("8th_order", "eighth_order") && return EighthOrderOcculter()
    throw(ArgumentError("Unknown occulter selector: $(occulter)"))
end

_coronagraph_plot_title(::GaussianOcculter) = "Gaussian spot"
_coronagraph_plot_title(::SolidOcculter) = "Solid spot"
_coronagraph_plot_title(::EighthOrderOcculter) = "8th order band limited spot"

function _plot_coronagraph_planes(after_occulter::AbstractMatrix, before_lyot::AbstractMatrix, occulter::AbstractOcculter)
    p1 = heatmap(
        after_occulter;
        aspect_ratio=:equal,
        color=:grays,
        colorbar=true,
        title="After Occulter",
    )
    p2 = heatmap(
        before_lyot;
        aspect_ratio=:equal,
        color=:grays,
        colorbar=true,
        title="Before Lyot Stop",
    )
    plt = plot(p1, p2; layout=(1, 2), size=(1200, 500), plot_title=_coronagraph_plot_title(occulter))
    display(plt)
    return plt
end

function _apply_occulter!(wfo::WaveFront, ::GaussianOcculter, occrad::Real, occrad_m::Real)
    r = prop_radius(wfo)
    h = sqrt(-0.5 * occrad_m^2 / log(1 - sqrt(0.5)))
    gauss_spot = @. 1 - exp(-0.5 * (r / h)^2)
    prop_multiply(wfo, gauss_spot)
    return wfo
end

function _apply_occulter!(wfo::WaveFront, ::SolidOcculter, occrad::Real, occrad_m::Real)
    prop_circular_obscuration(wfo, occrad_m)
    return wfo
end

function _apply_occulter!(wfo::WaveFront, ::EighthOrderOcculter, occrad::Real, occrad_m::Real)
    # `prop_8th_order_mask` interprets `occrad` in lambda/D units here.
    prop_8th_order_mask(wfo, occrad; circular=true)
    return wfo
end

function _apply_lyot_stop!(wfo::WaveFront, ::GaussianOcculter)
    prop_circular_aperture(wfo, 0.25; NORM=true)
    return wfo
end

function _apply_lyot_stop!(wfo::WaveFront, ::SolidOcculter)
    prop_circular_aperture(wfo, 0.84; NORM=true)
    return wfo
end

function _apply_lyot_stop!(wfo::WaveFront, ::EighthOrderOcculter)
    prop_circular_aperture(wfo, 0.50; NORM=true)
    return wfo
end

function coronagraph(wfo::WaveFront, f_lens::Real, occulter, diam::Real; PLOT::Bool=true, plot::Union{Nothing,Bool}=nothing)
    occulter = normalize_occulter(occulter)
    do_plot = plot === nothing ? PLOT : plot

    prop_lens(wfo, f_lens, "coronagraph imaging lens")
    prop_propagate(wfo, f_lens, "occulter")

    λ = prop_get_wavelength(wfo)
    occrad = 4.0
    occrad_rad = occrad * λ / diam
    dx_m = prop_get_sampling(wfo)
    dx_rad = prop_get_sampling_radians(wfo)
    occrad_m = occrad_rad * dx_m / dx_rad

    _apply_occulter!(wfo, occulter, occrad, occrad_m)

    after_occulter = do_plot ? sqrt.(prop_get_amplitude(wfo)) : nothing

    prop_propagate(wfo, f_lens, "pupil reimaging lens")
    prop_lens(wfo, f_lens, "pupil reimaging lens")
    prop_propagate(wfo, 2f_lens, "lyot stop")

    if do_plot
        before_lyot = prop_get_amplitude(wfo) .^ 0.2
        _plot_coronagraph_planes(after_occulter, before_lyot, occulter)
    end

    _apply_lyot_stop!(wfo, occulter)

    prop_propagate(wfo, f_lens, "reimaging lens")
    prop_lens(wfo, f_lens, "reimaging lens")
    prop_propagate(wfo, f_lens, "final focus")
    return wfo
end
