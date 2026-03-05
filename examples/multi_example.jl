using proper
using Plots

function multi_example(lambda_m::Real, n::Integer, passvalue=Dict("use_dm" => false, "dm" => zeros(48, 48)); kwargs...)
    diam = 0.048
    pupil_ratio = 0.25
    fl_lens = 0.48

    wfo = prop_begin(diam, lambda_m, n; beam_diam_fraction=pupil_ratio)
    prop_circular_aperture(wfo, diam / 2)
    prop_define_entrance(wfo)

    if get(passvalue, "use_dm", get(passvalue, :use_dm, false))
        dm = get(passvalue, "dm", get(passvalue, :dm, zeros(48, 48)))
        prop_dm(wfo, dm; mirror=false)
    end

    prop_lens(wfo, fl_lens)
    prop_propagate(wfo, fl_lens)
    return prop_end(wfo; noabs=true)
end

if abspath(PROGRAM_FILE) == @__FILE__
    stack, samplings = prop_run_multi(multi_example, 0.55, 256; PASSVALUE=[Dict("use_dm" => false), Dict("use_dm" => true, "dm" => zeros(48, 48))])
    println("multi_example: samplings = ", samplings)
    heatmap(log10.(abs.(stack[:, :, 1]) .+ eps()); aspect_ratio=:equal, title="multi_example [0]")
end
