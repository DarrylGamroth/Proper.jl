using Proper
using Plots
include(joinpath(@__DIR__, "_shared.jl"))
include(joinpath(@__DIR__, "run_coronagraph_dm.jl"))

function coronagraph_demo()
    n = 512
    lambda_um = 0.55

    no_errors, _ = prop_run(run_coronagraph_dm, lambda_um, n; PASSVALUE=Dict("use_errors" => false, "use_dm" => false, "occulter_type" => "8TH_ORDER"))
    with_errors, _ = prop_run(run_coronagraph_dm, lambda_um, n; PASSVALUE=Dict("use_errors" => true, "use_dm" => false, "occulter_type" => "8TH_ORDER"))
    with_dm, _ = prop_run(run_coronagraph_dm, lambda_um, n; PASSVALUE=Dict("use_errors" => true, "use_dm" => true, "occulter_type" => "8TH_ORDER"))

    nd = 256
    p1 = heatmap(center_crop(no_errors, nd) .^ 0.25; aspect_ratio=:equal, colorbar=false, title="No errors")
    p2 = heatmap(center_crop(with_errors, nd) .^ 0.25; aspect_ratio=:equal, colorbar=false, title="With errors")
    p3 = heatmap(center_crop(with_dm, nd) .^ 0.25; aspect_ratio=:equal, colorbar=false, title="DM corrected")
    display(plot(p1, p2, p3; layout=(1, 3), size=(1200, 400), plot_title="PSFs"))

    println("Maximum speckle flux / stellar flux:")
    println("  No wavefront errors = ", maximum(no_errors))
    println("  With wavefront errors = ", maximum(with_errors))
    println("  With DM correction = ", maximum(with_dm))

    return no_errors, with_errors, with_dm
end

if abspath(PROGRAM_FILE) == @__FILE__
    coronagraph_demo()
end
