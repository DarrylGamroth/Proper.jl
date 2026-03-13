using Test

@testset "Example ports smoke" begin
    exdir = joinpath(@__DIR__, "..", "examples")
    examples = (
        ("simple_prescription.jl", :simple_prescription),
        ("simple_telescope.jl", :simple_telescope),
        ("hubble_simple.jl", :hubble_simple),
        ("microscope.jl", :microscope),
        ("run_coronagraph.jl", :run_coronagraph),
        ("run_coronagraph_dm.jl", :run_coronagraph_dm),
        ("run_occulter.jl", :run_occulter),
        ("talbot.jl", :talbot),
        ("talbot_correct.jl", :talbot_correct),
        ("psdtest.jl", :psdtest),
        ("multi_example.jl", :multi_example),
    )

    for (file, sym) in examples
        mod = load_example_module(joinpath(exdir, file))
        @test isdefined(mod, sym)
    end
end
