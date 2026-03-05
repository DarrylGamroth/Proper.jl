using Test

@testset "Example ports smoke" begin
    exdir = joinpath(@__DIR__, "..", "examples")

    include(joinpath(exdir, "telescope.jl"))
    include(joinpath(exdir, "telescope_dm.jl"))
    include(joinpath(exdir, "coronagraph.jl"))
    include(joinpath(exdir, "run_occulter.jl"))
    include(joinpath(exdir, "run_coronagraph.jl"))
    include(joinpath(exdir, "run_coronagraph_dm.jl"))
    include(joinpath(exdir, "simple_prescription.jl"))
    include(joinpath(exdir, "simple_telescope.jl"))
    include(joinpath(exdir, "hubble_simple.jl"))
    include(joinpath(exdir, "microscope.jl"))
    include(joinpath(exdir, "psdtest.jl"))
    include(joinpath(exdir, "talbot.jl"))
    include(joinpath(exdir, "talbot_correct.jl"))
    include(joinpath(exdir, "multi_example.jl"))

    @test isdefined(Main, :simple_prescription)
    @test isdefined(Main, :simple_telescope)
    @test isdefined(Main, :hubble_simple)
    @test isdefined(Main, :microscope)
    @test isdefined(Main, :run_coronagraph)
    @test isdefined(Main, :run_coronagraph_dm)
    @test isdefined(Main, :run_occulter)
    @test isdefined(Main, :talbot)
    @test isdefined(Main, :talbot_correct)
    @test isdefined(Main, :psdtest)
    @test isdefined(Main, :multi_example)
end
