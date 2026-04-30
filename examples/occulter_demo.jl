using Proper
include(joinpath(@__DIR__, "run_occulter.jl"))

function occulter_demo()
    n = 512
    lambda_um = 0.55
    model = prepare_model(:run_occulter, run_occulter, lambda_um, n; pool_size=1)

    solid, _ = prop_run(model; occulter=:solid)
    gaussian, _ = prop_run(model; occulter=:gaussian)
    eighth_order, _ = prop_run(model; occulter=:eighth_order)

    return solid, gaussian, eighth_order
end

if abspath(PROGRAM_FILE) == @__FILE__
    solid, gaussian, eighth_order = occulter_demo()
    println("occulter_demo peaks: ", maximum(solid), ", ", maximum(gaussian), ", ", maximum(eighth_order))
end
