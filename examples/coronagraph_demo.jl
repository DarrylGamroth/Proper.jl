using Proper
using Plots
using Random
include(joinpath(@__DIR__, "_shared.jl"))
include(joinpath(@__DIR__, "run_coronagraph_dm.jl"))

function coronagraph_demo(
    ;
    n::Integer=512,
    lambda_um::Real=0.55,
    show_plot::Bool=true,
    map_file::AbstractString="telescope_obj.fits",
    rng::AbstractRNG=Random.default_rng(),
)
    n > 0 || throw(ArgumentError("n must be positive"))
    model = prepare_model(:run_coronagraph_dm, run_coronagraph_dm, lambda_um, n; pool_size=1)

    no_errors, _ = prop_run(
        model;
        use_errors=false,
        use_dm=false,
        occulter=:eighth_order,
        map_file,
        rng,
    )
    with_errors, _ = prop_run(
        model;
        use_errors=true,
        use_dm=false,
        occulter=:eighth_order,
        map_file,
        rng,
    )
    with_dm, _ = prop_run(
        model;
        use_errors=true,
        use_dm=true,
        occulter=:eighth_order,
        map_file,
        rng,
    )

    if show_plot
        nd = min(256, n)
        p1 = heatmap(center_crop(no_errors, nd) .^ 0.25; aspect_ratio=:equal, colorbar=false, title="No errors")
        p2 = heatmap(center_crop(with_errors, nd) .^ 0.25; aspect_ratio=:equal, colorbar=false, title="With errors")
        p3 = heatmap(center_crop(with_dm, nd) .^ 0.25; aspect_ratio=:equal, colorbar=false, title="DM corrected")
        display(plot(p1, p2, p3; layout=(1, 3), size=(1200, 400), plot_title="PSFs"))
    end

    println("Maximum speckle flux / stellar flux:")
    println("  No wavefront errors = ", maximum(no_errors))
    println("  With wavefront errors = ", maximum(with_errors))
    println("  With DM correction = ", maximum(with_dm))

    return no_errors, with_errors, with_dm
end

if abspath(PROGRAM_FILE) == @__FILE__
    coronagraph_demo()
end
