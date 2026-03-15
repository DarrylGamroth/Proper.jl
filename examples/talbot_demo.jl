using Proper
using Plots
using Statistics
include(joinpath(@__DIR__, "talbot.jl"))

function talbot_demo()
    diam = 0.1
    period = 0.04
    wavelength_microns = 0.5
    wavelength_m = wavelength_microns * 1e-6
    n = 128

    nseg = 9
    talbot_length = 2 * period^2 / wavelength_m
    delta_length = talbot_length / (nseg - 1)
    model = prepare_model(:talbot, talbot, wavelength_microns, n; pool_size=1)

    z = 0.0
    plots = Plot[]
    for _ in 1:nseg
        wavefront, _ = prop_run(model; PASSVALUE=Dict("diam" => diam, "period" => period, "dist" => z))
        line = wavefront[:, n ÷ 2 + 1]
        amp = abs.(line) .- mean(abs.(line))
        phase = angle.(line) .- mean(angle.(line))

        push!(plots, plot(amp; ylim=(-0.0015, 0.0015), title="Amplitude"))
        push!(plots, plot(phase; ylim=(-0.25, 0.25), title="Phase"))
        z += delta_length
    end

    display(plot(plots...; layout=(nseg, 2), size=(900, 1800)))
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    talbot_demo()
end
