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
        ("migration_dm_fits.jl", :migration_dm_fits_demo),
        ("wfirst_phaseb_reference.jl", :wfirst_phaseb_reference_demo),
    )

    for (file, sym) in examples
        mod = load_example_module(joinpath(exdir, file))
        @test isdefined(mod, sym)
    end
end

@testset "Migration DM/FITS example executes" begin
    mod = load_example_module(joinpath(@__DIR__, "..", "examples", "migration_dm_fits.jl"))
    psf, sampling = mod.migration_dm_fits_demo(; gridsize=32)
    @test size(psf) == (32, 32)
    @test sampling > 0
    @test maximum(psf) > 0
end

@testset "Multi example executes with direct DM map" begin
    mod = load_example_module(joinpath(@__DIR__, "..", "examples", "multi_example.jl"))
    field, sampling = mod.multi_example(0.55e-6, 32, Dict("use_dm" => true, "dm" => zeros(32, 32)))
    @test size(field) == (32, 32)
    @test sampling > 0
end

@testset "Example native keyword APIs remain compatible with PASSVALUE" begin
    exdir = joinpath(@__DIR__, "..", "examples")

    run_occulter_mod = load_example_module(joinpath(exdir, "run_occulter.jl"))
    native_solid, native_sampling = run_occulter_mod.run_occulter(0.55e-6, 32; occulter=:solid)
    compat_solid, compat_sampling = run_occulter_mod.run_occulter(0.55e-6, 32, Dict("occulter_type" => "SOLID"))
    @test native_solid ≈ compat_solid
    @test native_sampling == compat_sampling

    talbot_mod = load_example_module(joinpath(exdir, "talbot.jl"))
    native_field, native_talbot_sampling = talbot_mod.talbot(0.5e-6, 32; diam=0.1, period=0.04, dist=0.0)
    compat_field, compat_talbot_sampling = talbot_mod.talbot(0.5e-6, 32, Dict("diam" => 0.1, "period" => 0.04, "dist" => 0.0))
    @test native_field ≈ compat_field
    @test native_talbot_sampling == compat_talbot_sampling
end
