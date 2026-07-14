using Proper
using Random
using Test

function load_example_module(path::AbstractString; name::Symbol=gensym(Symbol(splitext(basename(path))[1])))
    mod = Module(name)
    Core.eval(mod, :(include(path::AbstractString) = Base.include($mod, path)))
    Base.include(mod, path)
    return mod
end

const EXAMPLE_DIR = @__DIR__
const UPSTREAM_EXAMPLE_EXPORTS = (
    ("__init__.jl", :ProperExamples),
    ("coronagraph.jl", :coronagraph),
    ("coronagraph_demo.jl", :coronagraph_demo),
    ("example_system.jl", :example_system),
    ("hubble_simple.jl", :hubble_simple),
    ("microscope.jl", :microscope),
    ("multi_example.jl", :multi_example),
    ("occulter_demo.jl", :occulter_demo),
    ("psdtest.jl", :psdtest),
    ("run_coronagraph.jl", :run_coronagraph),
    ("run_coronagraph_dm.jl", :run_coronagraph_dm),
    ("run_example.jl", :run_example),
    ("run_occulter.jl", :run_occulter),
    ("simple_prescription.jl", :simple_prescription),
    ("simple_telescope.jl", :simple_telescope),
    ("talbot.jl", :talbot),
    ("talbot_correct.jl", :talbot_correct),
    ("talbot_correct_demo.jl", :talbot_correct_demo),
    ("talbot_demo.jl", :talbot_demo),
    ("telescope.jl", :telescope),
    ("telescope_dm.jl", :telescope_dm),
    ("testmulti1.jl", :testmulti1),
    ("testmulti2.jl", :testmulti2),
)

@testset "Upstream example file and entry-point coverage" begin
    @test length(UPSTREAM_EXAMPLE_EXPORTS) == 23
    for (file, entrypoint) in UPSTREAM_EXAMPLE_EXPORTS
        path = joinpath(EXAMPLE_DIR, file)
        @test isfile(path)
        mod = load_example_module(path)
        @test isdefined(mod, entrypoint)
    end

    for file in filter(
        path -> endswith(path, ".jl") && basename(path) != "runtests.jl",
        readdir(EXAMPLE_DIR; join=true),
    )
        source = read(file, String)
        @test !occursin("Any[]", source)
        @test !occursin("Vector{Any}", source)
    end
end

@testset "Julia-specific examples load" begin
    for (file, entrypoint) in (
        ("migration_dm_fits.jl", :migration_dm_fits_demo),
        ("wfirst_phaseb_reference.jl", :wfirst_phaseb_reference_demo),
    )
        mod = load_example_module(joinpath(EXAMPLE_DIR, file))
        @test isdefined(mod, entrypoint)
    end
end

@testset "Example runners execute without display or repository writes" begin
    migration_mod = load_example_module(joinpath(EXAMPLE_DIR, "migration_dm_fits.jl"))
    migration_psf, migration_sampling = migration_mod.migration_dm_fits_demo(; gridsize=32)
    @test size(migration_psf) == (32, 32)
    @test all(isfinite, migration_psf)
    @test migration_sampling > 0

    multi_mod = load_example_module(joinpath(EXAMPLE_DIR, "multi_example.jl"))
    no_dm, no_dm_sampling = multi_mod.multi_example(0.55e-6, 32; use_dm=false)
    zero_dm, zero_dm_sampling = multi_mod.multi_example(
        0.55e-6,
        32;
        use_dm=true,
        dm=zeros(48, 48),
    )
    @test zero_dm == no_dm
    @test zero_dm_sampling == no_dm_sampling

    testmulti1_mod = load_example_module(joinpath(EXAMPLE_DIR, "testmulti1.jl"))
    broadband_psf = testmulti1_mod.testmulti1(; gridsize=32, npsf=32, nlambda=2)
    @test size(broadband_psf) == (32, 32)
    @test all(isfinite, broadband_psf)
    @test maximum(broadband_psf) > 0

    testmulti2_mod = load_example_module(joinpath(EXAMPLE_DIR, "testmulti2.jl"))
    ripple_fields, ripple_samplings = testmulti2_mod.testmulti2(; gridsize=32, npatterns=2)
    @test size(ripple_fields) == (32, 32, 2)
    @test all(isfinite, ripple_fields)
    @test length(ripple_samplings) == 2
    @test all(>(0), ripple_samplings)

    run_example_mod = load_example_module(joinpath(EXAMPLE_DIR, "run_example.jl"))
    state_psf, state_sampling = run_example_mod.run_example(0.5e-6, 32; iterations=2)
    @test size(state_psf) == (32, 32)
    @test all(isfinite, state_psf)
    @test state_sampling > 0

    occulter_mod = load_example_module(joinpath(EXAMPLE_DIR, "occulter_demo.jl"))
    occulter_outputs = occulter_mod.occulter_demo(; n=32)
    @test all(output -> size(output) == (32, 32), occulter_outputs)
    @test all(output -> all(isfinite, output), occulter_outputs)

    talbot_mod = load_example_module(joinpath(EXAMPLE_DIR, "talbot_demo.jl"))
    @test isnothing(talbot_mod.talbot_demo(; n=32, nseg=3, show_plot=false))

    talbot_correct_mod = load_example_module(joinpath(EXAMPLE_DIR, "talbot_correct_demo.jl"))
    @test isnothing(talbot_correct_mod.talbot_correct_demo(; n=32, nseg=3, show_plot=false))

    coronagraph_mod = load_example_module(joinpath(EXAMPLE_DIR, "coronagraph_demo.jl"))
    mktempdir() do dir
        outputs = coronagraph_mod.coronagraph_demo(
            ;
            n=64,
            show_plot=false,
            map_file=joinpath(dir, "telescope_obj.fits"),
            rng=MersenneTwister(0x5cc),
        )
        @test all(output -> size(output) == (64, 64), outputs)
        @test all(output -> all(isfinite, output), outputs)
    end
end
