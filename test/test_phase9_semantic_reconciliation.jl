using Test
using Random

@testset "Phase 9 semantic reconciliation hotspots" begin
    @testset "prop_resamplemap uses independent yshift" begin
        wf = prop_begin(1.0, 500e-9, 16)
        dmap = reshape(collect(1.0:256.0), 16, 16)
        pix = wf.sampling_m

        m0 = prop_resamplemap(wf, dmap, pix, 8.0, 8.0, 0.0, 0.0)
        my = prop_resamplemap(wf, dmap, pix, 8.0, 8.0, 0.0, 2pix)
        mx = prop_resamplemap(wf, dmap, pix, 8.0, 8.0, 2pix, 0.0)

        @test !isapprox(my, m0; atol=0, rtol=0)
        @test !isapprox(mx, m0; atol=0, rtol=0)
        @test !isapprox(mx, my; atol=0, rtol=0)
    end

    @testset "prop_end extract semantics are integer-safe" begin
        wf = prop_begin(1.0, 500e-9, 15)
        out5, _ = prop_end(wf; extract=5)
        @test size(out5) == (5, 5)

        wf2 = prop_begin(1.0, 500e-9, 16)
        out6, _ = prop_end(wf2; extract=6)
        @test size(out6) == (6, 6)
    end

    @testset "prop_state restores full wavefront state" begin
        wf = prop_begin(1.0, 500e-9, 16)
        wf.field .*= cis.(rand(MersenneTwister(7), 16, 16))
        wf.wavelength_m = 632.8e-9
        wf.sampling_m = 2.3e-4
        wf.z_m = 0.456
        wf.beam_diameter_m = 0.321
        wf.z_w0_m = 0.654
        wf.w0_m = 0.123
        wf.z_rayleigh_m = 1.234
        wf.current_fratio = 42.0
        wf.reference_surface = Proper.SPHERICAL
        wf.beam_type_old = Proper.OUTSIDE
        wf.propagator_type = Proper.OUTSIDE_TO_OUTSIDE
        wf.rayleigh_factor = 1.5

        mktempdir() do d
            state = joinpath(d, "wf.state")
            prop_savestate(wf, state)

            wf2 = prop_begin(1.0, 500e-9, 16)
            prop_state(wf2, state)

            @test all(isapprox.(wf2.field, wf.field))
            @test wf2.wavelength_m == wf.wavelength_m
            @test wf2.sampling_m == wf.sampling_m
            @test wf2.z_m == wf.z_m
            @test wf2.beam_diameter_m == wf.beam_diameter_m
            @test wf2.z_w0_m == wf.z_w0_m
            @test wf2.w0_m == wf.w0_m
            @test wf2.z_rayleigh_m == wf.z_rayleigh_m
            @test wf2.current_fratio == wf.current_fratio
            @test wf2.reference_surface == wf.reference_surface
            @test wf2.beam_type_old == wf.beam_type_old
            @test wf2.propagator_type == wf.propagator_type
            @test wf2.rayleigh_factor == wf.rayleigh_factor
        end
    end

    @testset "prop_psd_errormap FILE reuse semantics" begin
        wf = prop_begin(1.0, 500e-9, 32)
        mktempdir() do d
            f = joinpath(d, "psd_map.fits")
            m1 = prop_psd_errormap(wf, 1e-18, 10.0, 3.0; FILE=f, no_apply=true, rng=MersenneTwister(123))
            @test isfile(f)

            # Different parameters should still reuse existing FILE map.
            m2 = prop_psd_errormap(wf, 9e-18, 1.0, 1.5; FILE=f, no_apply=true, rng=MersenneTwister(999))
            @test size(m2) == size(m1)
            @test maximum(abs.(m2 .- m1)) < 1e-6
        end
    end

    @testset "DM correction behavior matches upstream trend" begin
        exdir = joinpath(@__DIR__, "..", "examples")
        include(joinpath(exdir, "telescope.jl"))
        include(joinpath(exdir, "telescope_dm.jl"))
        include(joinpath(exdir, "coronagraph.jl"))
        include(joinpath(exdir, "run_coronagraph_dm.jl"))

        mktempdir() do d
            cp(joinpath(exdir, "telescope_obj.fits"), joinpath(d, "telescope_obj.fits"); force=true)
            cd(d) do
                psf_err, _ = run_coronagraph_dm(0.55e-6, 256, Dict("use_errors" => true, "use_dm" => false, "occulter_type" => "GAUSSIAN"))
                psf_dm, _ = run_coronagraph_dm(0.55e-6, 256, Dict("use_errors" => true, "use_dm" => true, "occulter_type" => "GAUSSIAN"))
                @test maximum(psf_dm) / maximum(psf_err) < 0.01
            end
        end
    end
end
