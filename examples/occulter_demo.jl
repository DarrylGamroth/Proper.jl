using Proper
include(joinpath(@__DIR__, "run_occulter.jl"))

function occulter_demo()
    n = 512
    lambda_um = 0.55

    solid, _ = prop_run(run_occulter, lambda_um, n; PASSVALUE=Dict("occulter_type" => "SOLID"))
    gaussian, _ = prop_run(run_occulter, lambda_um, n; PASSVALUE=Dict("occulter_type" => "GAUSSIAN"))
    eighth_order, _ = prop_run(run_occulter, lambda_um, n; PASSVALUE=Dict("occulter_type" => "8TH_ORDER"))

    return solid, gaussian, eighth_order
end

if abspath(PROGRAM_FILE) == @__FILE__
    solid, gaussian, eighth_order = occulter_demo()
    println("occulter_demo peaks: ", maximum(solid), ", ", maximum(gaussian), ", ", maximum(eighth_order))
end
