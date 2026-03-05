using proper
using Plots
using Statistics
include(joinpath(@__DIR__, "_shared.jl"))

function talbot(wavelength::Real, gridsize::Integer, passvalue=Dict("period" => 0.0, "diam" => 0.0, "dist" => 0.0))
    period = get(passvalue, "period", get(passvalue, :period, 0.0))
    diam = get(passvalue, "diam", get(passvalue, :diam, 0.0))
    dist = get(passvalue, "dist", get(passvalue, :dist, 0.0))

    wfo = prop_begin(diam, wavelength, gridsize)
    m = 0.2
    x = (collect(0:(gridsize - 1)) .- gridsize / 2) .* prop_get_sampling(wfo)
    grating_1d = @. 0.5 * (1 + m * cos(2pi * x / period))
    grating = grating_1d .* ones(Float64, 1, gridsize)

    prop_multiply(wfo, grating)
    prop_define_entrance(wfo)
    prop_propagate(wfo, dist; TO_PLANE=true)
    return prop_end(wfo; noabs=true)
end

if abspath(PROGRAM_FILE) == @__FILE__
    wavefront, sampling = talbot(0.5e-6, 128, Dict("period" => 0.04, "diam" => 0.1, "dist" => 0.0))
    println("talbot: sampling = ", sampling)
    plot(abs.(wavefront[:, size(wavefront, 2) ÷ 2 + 1]) .- mean(abs.(wavefront[:, size(wavefront, 2) ÷ 2 + 1])); title="talbot amplitude")
end
