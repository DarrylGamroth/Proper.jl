using Proper
using Plots
using Statistics
include(joinpath(@__DIR__, "talbot.jl"))

function talbot_demo(
    ;
    diam::Real=0.1,
    period::Real=0.04,
    wavelength_microns::Real=0.5,
    n::Integer=128,
    nseg::Integer=9,
    show_plot::Bool=true,
)
    n > 0 || throw(ArgumentError("n must be positive"))
    nseg > 1 || throw(ArgumentError("nseg must be at least two"))
    wavelength_m = wavelength_microns * 1e-6

    talbot_length = 2 * period^2 / wavelength_m
    delta_length = talbot_length / (nseg - 1)
    model = prepare_model(:talbot, talbot, wavelength_microns, n; pool_size=1)

    profiles = map(range(0.0; step=delta_length, length=nseg)) do z
        wavefront, _ = prop_run(model; diam=diam, period=period, dist=z)
        line = wavefront[:, n ÷ 2 + 1]
        amp = abs.(line) .- mean(abs.(line))
        phase = angle.(line) .- mean(angle.(line))
        return (; amplitude=amp, phase)
    end

    if show_plot
        plots = [
            plot(series; ylim=limits, title)
            for profile in profiles
            for (series, limits, title) in (
                (profile.amplitude, (-0.0015, 0.0015), "Amplitude"),
                (profile.phase, (-0.25, 0.25), "Phase"),
            )
        ]
        display(plot(plots...; layout=(nseg, 2), size=(900, 1800)))
    end
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    talbot_demo()
end
