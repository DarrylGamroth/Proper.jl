using Test
using Random
using Statistics

@testset "Phase 2 core kernels" begin
    wf = prop_begin(1.0, 500e-9, 32)
    orig = copy(wf.field)

    prop_lens(wf, 10.0)
    @test wf.field != orig
    @test isfinite(prop_get_fratio(wf))

    wf2 = prop_begin(1.0, 500e-9, 32)
    prop_circular_aperture(wf2, 0.2)
    @test any(iszero, abs.(wf2.field))

    wf3 = prop_begin(1.0, 500e-9, 32)
    prop_propagate(wf3, 0.25)
    @test isapprox(prop_get_z(wf3), 0.25; atol=1e-12)
    @test wf3.reference_surface in (:PLANAR, :SPHERI)

    a = reshape(collect(1.0:16.0), 4, 4)
    @test size(prop_magnify(a, 2.0)) == (8, 8)
    @test size(prop_rotate(a, 15.0)) == size(a)

    wf4 = prop_begin(1.0, 500e-9, 32)
    @test (@inferred prop_select_propagator(wf4, 0.1)) isa Float64
    @test (@inferred prop_qphase(wf4, 1.2)) === wf4
    @test (@inferred prop_ptp(wf4, 0.01)) === wf4
end

@testset "Phase 2 run entrypoints" begin
    dummy(λm, n; kwargs...) = prop_begin(1.0, λm, n)

    out, s = prop_run(dummy, 0.55, 16)
    @test size(out) == (16, 16)
    @test s > 0

    passvals = [1, 2, 3]
    stack, samplings = prop_run_multi(dummy, 0.55, 16; PASSVALUE=passvals)
    @test size(stack) == (16, 16, 3)
    @test length(samplings) == 3
end

@testset "Phase 3 FITS/map basics" begin
    wf = prop_begin(1.0, 500e-9, 16)
    dmap = randn(16, 16)

    mktempdir() do d
        f = joinpath(d, "map.fits")
        prop_writemap(dmap, f; SAMPLING=wf.sampling_m)
        m2 = prop_readmap(wf, f; SAMPLING=wf.sampling_m)
        @test size(m2) == size(wf.field)

        f2 = joinpath(d, "img.fits")
        prop_fits_write(f2, dmap)
        r2 = prop_fits_read(f2)
        @test size(r2) == size(dmap)
    end
end

@testset "Phase 4 state and compatibility helpers" begin
    wf = prop_begin(1.0, 500e-9, 16)
    wf.field .*= 2
    wf.z_m = 0.123

    mktempdir() do d
        sfile = joinpath(d, "wf.state")
        @test !prop_is_statesaved(sfile)
        prop_savestate(wf, sfile)
        @test prop_is_statesaved(sfile)

        wf2 = prop_begin(1.0, 500e-9, 16)
        prop_state(wf2, sfile)
        @test isapprox(prop_get_z(wf2), 0.123; atol=1e-12)
        @test all(isapprox.(wf2.field, wf.field))

        out, _ = prop_end_savestate(wf2, joinpath(d, "wf2.state"))
        @test size(out) == (16, 16)
    end

    @test prop_use_fftw() == true
    @test prop_use_ffti() == false
    @test prop_set_antialiasing(3) == 3.0
end

@testset "Phase 5 placeholder API coverage" begin
    wf = prop_begin(1.0, 500e-9, 32)
    dm = randn(32, 32) .* 1e-9
    @test prop_dm(wf, dm) === wf
    @test prop_dm(wf, dm; mirror=true) === wf
    @test_throws MethodError prop_dm(wf, dm; MIRROR=1)

    @test prop_sinc(0.0) == 1.0
    @test length(prop_noll_zernikes(5)) == 5
    @test size(prop_zernikes(wf, 3)) == (32, 32, 3)
    @test length(prop_fit_zernikes(wf, 3)) == 3

    pix = prop_pixellate(rand(32, 32), 2)
    @test size(pix) == (16, 16)

    m1 = prop_polygon(wf, 6, 0.2)
    m2 = prop_irregular_polygon(wf, [-0.1, 0.1, 0.1, -0.1], [-0.1, -0.1, 0.1, 0.1])
    m3 = prop_rounded_rectangle(wf, 0.05, 0.2, 0.3)
    @test size(m1) == size(wf.field)
    @test size(m2) == size(wf.field)
    @test size(m3) == size(wf.field)

    @test prop_fftw() == true
    @test prop_ffti() == false
    @test prop_compile_c() === nothing
    @test prop_dftidefs()[:DFTI_FORWARD_SCALE] == 1.0
    @test libcconv(rand(8, 8), 4.2, 3.7) isa Real
    @test libcconvthread(rand(8, 8), 4.2, 3.7) isa Real
    @test size(libszoom(rand(8, 8), 2.0), 1) == 16
    @test size(prop_szoom(rand(8, 8), 2.0), 1) == 16
end

@testset "Phase 8 geometry and zernike parity groundwork" begin
    wf = prop_begin(1.0, 500e-9, 64)

    # Noll ordering sanity
    desc = prop_noll_zernikes(6)
    @test desc[1] == (n=0, m=0, trig=:none)
    @test desc[2].n == 1 && desc[2].m == 1
    @test desc[4] == (n=2, m=0, trig=:none)

    # Zernike map generation and application
    zmap = prop_zernikes(wf, [4, 5], [1e-9, 2e-9]; no_apply=true)
    @test size(zmap) == size(wf.field)

    coeff, fitmap = prop_fit_zernikes(zmap, ones(size(zmap)...), 32.0, 6; fit=true)
    @test length(coeff) == 6
    @test size(fitmap) == size(zmap)

    # Polygon/irregular polygon/rounded rectangle masks should be non-trivial.
    mpoly = prop_polygon(wf, 6, 0.2)
    mirr = prop_irregular_polygon(wf, [-0.2, 0.2, 0.2, -0.2], [-0.2, -0.2, 0.2, 0.2])
    mround = prop_rounded_rectangle(wf, 0.05, 0.3, 0.2)
    @test 0.0 < mean(mpoly) < 1.0
    @test 0.0 < mean(mirr) < 1.0
    @test 0.0 < mean(mround) < 1.0

    # PSD map generation should produce finite-valued map.
    dmap = prop_psd_errormap(wf, 1e-18, 10.0, 3.0; no_apply=true, rng=MersenneTwister(1))
    @test size(dmap) == size(wf.field)
    @test all(isfinite, dmap)
end

@testset "Phase 8 segmentation and interpolation parity coverage" begin
    wf = prop_begin(1.0, 500e-9, 64)
    m = prop_8th_order_mask(wf, 3.0; circular=true)
    @test size(m) == size(wf.field)
    @test 0.0 <= minimum(m) <= maximum(m) <= 1.0

    hz = prop_hex_zernikes(1:3, Float32[1e-9, -2e-9, 1e-9], 64, Float32(prop_get_sampling(wf)), 0.03f0)
    @test size(hz) == (64, 64)
    @test eltype(hz) == Float32

    nhex = 1 * (1 + 1) * 3 + 1
    zvals = zeros(Float64, 22, nhex)
    zvals[4, :] .= 1e-9
    ap, ph = prop_hex_wavefront(wf, 1, 0.03, 0.032, zvals; no_apply=true)
    @test size(ap) == size(wf.field)
    @test size(ph) == size(wf.field)
    @test any(!iszero, ap)

    a = reshape(collect(1.0:16.0), 4, 4)
    @test prop_cubic_conv(a, 2.0, 2.0) ≈ a[2, 2] atol=1e-12
    out = prop_cubic_conv(a, [1.0, 2.0, 3.0], [1.0, 2.0]; grid=true)
    @test size(out) == (2, 3)
    @test size(libszoom(a, 2.0), 1) == 8
end
