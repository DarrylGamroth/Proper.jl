using Test
using Random

@testset "R5 expanded inference/allocation gates" begin
    rng = MersenneTwister(2026)

    @testset "Propagation core kernels" begin
        wf_select = prop_begin(1.0, 500e-9, 64)
        @test (@inferred prop_select_propagator(wf_select, 0.05)) isa Float64
        @test (@allocated prop_select_propagator(wf_select, 0.01)) == 0

        wf_q = prop_begin(1.0, 500e-9, 64)
        @test (@inferred prop_qphase(wf_q, 0.1)) === wf_q
        prop_qphase(wf_q, 0.1) # warmup
        @test (@allocated prop_qphase(wf_q, 0.1)) < 300_000

        wf_ptp = prop_begin(1.0, 500e-9, 64)
        ctx_ptp = RunContext(typeof(wf_ptp.field))
        @test (@inferred prop_ptp(wf_ptp, 0.01, ctx_ptp)) === wf_ptp
        prop_ptp(wf_ptp, 0.01, ctx_ptp) # warmup
        @test (@allocated prop_ptp(wf_ptp, 0.01, ctx_ptp)) < 6_000_000

        wf_lens = prop_begin(1.0, 500e-9, 64)
        @test (@inferred prop_lens(wf_lens, 10.0)) === wf_lens
        prop_lens(wf_lens, 10.0) # warmup
        @test (@allocated prop_lens(wf_lens, 10.0)) < 35_000_000

        wf_prop = prop_begin(1.0, 500e-9, 64)
        ctx_prop = RunContext(typeof(wf_prop.field))
        @test (@inferred prop_propagate(wf_prop, 0.01, ctx_prop)) === wf_prop
        prop_propagate(wf_prop, 0.01, ctx_prop) # warmup
        @test (@allocated prop_propagate(wf_prop, 0.01, ctx_prop)) < 2_000_000
    end

    @testset "Map/PSD hotspots" begin
        wf = prop_begin(1.0, 500e-9, 64)
        ctx = RunContext(typeof(wf.field))
        dmap = rand(rng, Float32, 64, 64)

        opts = Proper.ResampleMapOptions(wf, wf.sampling_m, 32.0, 32.0)
        out = zeros(Float64, size(wf.field)...)
        @test (@inferred Proper.prop_resamplemap!(out, wf, dmap, opts, ctx)) === out
        Proper.prop_resamplemap!(out, wf, dmap, opts, ctx) # warmup
        @test (@allocated Proper.prop_resamplemap!(out, wf, dmap, opts, ctx)) == 0

        img = rand(rng, Float32, 64, 64)
        rot_out = similar(img)
        @test (@inferred prop_rotate!(rot_out, img, 9.0, ctx)) === rot_out
        prop_rotate!(rot_out, img, 9.0, ctx) # warmup
        @test (@allocated prop_rotate!(rot_out, img, 9.0, ctx)) < 250_000

        mag_out = similar(img, 64, 64)
        @test (@inferred prop_magnify!(mag_out, img, 1.1, ctx; QUICK=true)) === mag_out
        prop_magnify!(mag_out, img, 1.1, ctx; QUICK=true) # warmup
        @test (@allocated prop_magnify!(mag_out, img, 1.1, ctx; QUICK=true)) < 20_000

        wf_psd = prop_begin(0.212, 500e-9, 64)
        m = prop_psd_errormap(wf_psd, 3.29e-23, 212.26, 7.8; no_apply=true, rng=MersenneTwister(7))
        @test m isa AbstractMatrix
        @test size(m) == size(wf_psd.field)
        prop_psd_errormap(wf_psd, 3.29e-23, 212.26, 7.8; no_apply=true, rng=MersenneTwister(7)) # warmup
        @test (@allocated prop_psd_errormap(wf_psd, 3.29e-23, 212.26, 7.8; no_apply=true, rng=MersenneTwister(7))) < 500_000_000
    end

    @testset "Geometry mask hotspots" begin
        wf = prop_begin(1.0, 500e-9, 64)
        RT = real(eltype(wf.field))
        n = size(wf.field, 1)

        ellipse_out = zeros(RT, n, n)
        @test (@inferred Proper.prop_ellipse!(ellipse_out, wf, 0.35, 0.2, 0.0, 0.0; NORM=true)) === ellipse_out
        Proper.prop_ellipse!(ellipse_out, wf, 0.35, 0.2, 0.0, 0.0; NORM=true) # warmup
        @test (@allocated Proper.prop_ellipse!(ellipse_out, wf, 0.35, 0.2, 0.0, 0.0; NORM=true)) < 50_000

        rect_out = zeros(RT, n, n)
        @test (@inferred Proper.prop_rectangle!(rect_out, wf, 0.4, 0.2, 0.0, 0.0; NORM=true)) === rect_out
        Proper.prop_rectangle!(rect_out, wf, 0.4, 0.2, 0.0, 0.0; NORM=true) # warmup
        @test (@allocated Proper.prop_rectangle!(rect_out, wf, 0.4, 0.2, 0.0, 0.0; NORM=true)) < 50_000

        poly_out = zeros(RT, n, n)
        @test (@inferred Proper.prop_polygon!(poly_out, wf, 6, 0.33, 0.0, 0.0; NORM=true)) === poly_out
        Proper.prop_polygon!(poly_out, wf, 6, 0.33, 0.0, 0.0; NORM=true) # warmup
        @test (@allocated Proper.prop_polygon!(poly_out, wf, 6, 0.33, 0.0, 0.0; NORM=true)) < 100_000

        xverts = RT[-0.20, 0.12, 0.28, -0.08]
        yverts = RT[-0.18, -0.22, 0.19, 0.25]
        ipoly_out = zeros(RT, n, n)
        @test (@inferred Proper.prop_irregular_polygon!(ipoly_out, wf, xverts, yverts; NORM=true)) === ipoly_out
        Proper.prop_irregular_polygon!(ipoly_out, wf, xverts, yverts; NORM=true) # warmup
        @test (@allocated Proper.prop_irregular_polygon!(ipoly_out, wf, xverts, yverts; NORM=true)) < 100_000
    end
end
