using Test
using Random
using Statistics
using Base.Threads
using FFTW

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
    @test wf3.reference_surface in (Proper.PLANAR, Proper.SPHERICAL)

    a = reshape(collect(1.0:16.0), 4, 4)
    @test size(prop_magnify(a, 2.0)) == (8, 8)
    @test size(prop_rotate(a, 15.0)) == size(a)

    wf4 = prop_begin(1.0, 500e-9, 32)
    @test (@inferred prop_select_propagator(wf4, 0.1)) isa Float64
    @test (@inferred prop_qphase(wf4, 1.2)) === wf4
    @test (@inferred prop_ptp(wf4, 0.01)) === wf4
end

@testset "Centered core helpers" begin
    a = reshape(collect(1.0:25.0), 5, 5)
    out_crop = Matrix{Float64}(undef, 4, 4)
    Proper.copy_centered!(out_crop, a)
    @test out_crop == [1.0 6.0 11.0 16.0; 2.0 7.0 12.0 17.0; 3.0 8.0 13.0 18.0; 4.0 9.0 14.0 19.0]

    out_pad = Matrix{Float64}(undef, 7, 7)
    Proper.copy_centered!(out_pad, a)
    ref_pad = zeros(7, 7)
    ref_pad[2:6, 2:6] .= a
    @test out_pad == ref_pad

    c = reshape(ComplexF64.(1:16), 4, 4)
    shifted = similar(c)
    Proper.half_shift_copy!(shifted, c)
    @test shifted == prop_shift_center(c)

    odd = reshape(ComplexF64.(1:25), 5, 5)
    odd_shifted = similar(odd)
    Proper.shift_copy!(odd_shifted, odd, fld(size(odd, 1), 2), fld(size(odd, 2), 2))
    @test odd_shifted == prop_shift_center(odd)
    Proper.shift_copy!(odd_shifted, odd, cld(size(odd, 1), 2), cld(size(odd, 2), 2))
    @test odd_shifted == prop_shift_center(odd; inverse=true)

    rev = similar(c)
    Proper.reverse_shift1!(rev, c)
    @test rev == circshift(reverse(reverse(c; dims=1); dims=2), (1, 1))

    cache = Proper.CenteredFFTCache(Float64, 4; flags=FFTW.ESTIMATE)
    fwd = copy(c)
    Proper.centered_fft!(fwd, cache, -1)
    @test fwd ≈ circshift(fft(circshift(c, (-2, -2))) ./ length(c), (2, 2))

    invf = copy(c)
    Proper.centered_fft!(invf, cache, +1)
    @test invf ≈ circshift(ifft(circshift(c, (-2, -2))) .* length(c), (2, 2))
end

@testset "Odd-grid shift application semantics" begin
    map = reshape(collect(1.0:25.0), 5, 5)

    wf_mul = prop_begin(1.0, 500e-9, 5)
    fill!(wf_mul.field, 1 + 0im)
    prop_multiply(wf_mul, map)
    @test wf_mul.field == ComplexF64.(FFTW.ifftshift(map))

    wf_div = prop_begin(1.0, 500e-9, 5)
    fill!(wf_div.field, 1 + 0im)
    prop_divide(wf_div, map)
    @test wf_div.field ≈ ComplexF64.(inv.(FFTW.ifftshift(map)))

    phase_map = zeros(5, 5)
    phase_map[3, 3] = 1e-9
    wf_phase = prop_begin(1.0, 500e-9, 5)
    fill!(wf_phase.field, 1 + 0im)
    prop_add_phase(wf_phase, phase_map)
    @test wf_phase.field ≈ cis.((2pi / wf_phase.wavelength_m) .* FFTW.ifftshift(phase_map))
end

@testset "Phase 2 run entrypoints" begin
    dummy(λm, n; kwargs...) = prop_begin(1.0, λm, n)

    out, s = prop_run(dummy, 0.55, 16)
    @test size(out) == (16, 16)
    @test s > 0

    prepared = prepare_prescription(dummy, 0.55, 16)
    out_prepared, s_prepared = prop_run(prepared)
    @test size(out_prepared) == (16, 16)
    @test s_prepared > 0

    passvals = [1, 2, 3]
    stack, samplings = prop_run_multi(dummy, 0.55, 16; PASSVALUE=passvals)
    @test size(stack) == (16, 16, 3)
    @test length(samplings) == 3

    prepared_multi = prepare_prescription(dummy, 0.55, 16)
    stack_prepared, samplings_prepared = prop_run_multi(prepared_multi; PASSVALUE=passvals)
    @test size(stack_prepared) == (16, 16, 3)
    @test length(samplings_prepared) == 3

    prepared_batch = prepare_prescription_batch(prepared_multi; pool_size=2)
    stack_batch, samplings_batch = prop_run_multi(prepared_batch; PASSVALUE=passvals)
    @test size(stack_batch) == (16, 16, 3)
    @test length(samplings_batch) == 3

    prepared_model = prepare_model(dummy, 0.55, 16; pool_size=2)
    out_model, s_model = prop_run(prepared_model)
    @test size(out_model) == (16, 16)
    @test s_model > 0
    stack_model, samplings_model = prop_run_multi(prepared_model; PASSVALUE=passvals)
    @test size(stack_model) == (16, 16, 3)
    @test length(samplings_model) == 3

    dummy_pass(λm, n, pass; kwargs...) = prop_begin(1.0 + pass, λm, n)
    prepared_pass = prepare_prescription(dummy_pass, 0.55, 16; PASSVALUE=2)
    out_pass, s_pass = prop_run(prepared_pass)
    @test size(out_pass) == (16, 16)
    @test s_pass > 0

    let ctx = RunContext(Matrix{Float32})
        function dummy_scoped(λm, n; kwargs...)
            wf = prop_begin(1.0, λm, n)
            @test eltype(wf.field) == ComplexF32
            @test wf.workspace === ctx.workspace
            return wf
        end

        out_ctx, s_ctx = prop_run(dummy_scoped, 0.55f0, 16; context=ctx)
        @test size(out_ctx) == (16, 16)
        @test eltype(out_ctx) == Float32
        @test s_ctx > 0

        prepared_ctx = prepare_prescription(dummy_scoped, 0.55f0, 16; context=ctx)
        out_prepared_ctx, s_prepared_ctx = prop_run(prepared_ctx)
        @test size(out_prepared_ctx) == (16, 16)
        @test eltype(out_prepared_ctx) == Float32
        @test s_prepared_ctx > 0

        prepared_ctxs = Proper.prepared_contexts(prepared_ctx, 3)
        @test length(prepared_ctxs) == 3
        @test all(c -> c isa RunContext, prepared_ctxs)
        @test all(c -> c !== ctx, prepared_ctxs)
        @test all(c -> c.workspace !== ctx.workspace, prepared_ctxs)
        @test all(c -> eltype(Proper.field_backend_template(c.workspace)) == ComplexF32, prepared_ctxs)

        prepared_batch_ctx = prepare_prescription_batch(prepared_ctx; pool_size=2)
        @test length(prepared_batch_ctx.contexts) == 2
        Proper.ensure_prepared_batch_contexts!(prepared_batch_ctx, 4)
        @test length(prepared_batch_ctx.contexts) == 4
        @test all(c -> c.workspace !== ctx.workspace, prepared_batch_ctx.contexts)

        prepared_model_ctx = prepare_model(prepared_ctx; name=:dummy_scoped, pool_size=2)
        @test prepared_model_ctx.name == :dummy_scoped
        @test Proper.prepared_assets(prepared_model_ctx) === nothing
        @test Proper.prepared_batch(prepared_model_ctx) === prepared_model_ctx.batch
        @test Proper.prepared_prescription(prepared_model_ctx) === prepared_model_ctx.prepared

        dummy_prepared_parallel(λm, n, pass; kwargs...) = begin
            active = Proper.active_run_context()
            @test active isa RunContext
            if haskey(kwargs, :asset_slot)
                @test kwargs[:asset_slot] == pass
                @test kwargs[:asset_model] == "asset_model"
            end
            psf = fill(Float32(pass), n, n)
            psf[1, 1] = Float32(objectid(active.workspace) % 1024)
            return psf, 1.0f0
        end

        asset_counter = Ref(0)
        asset_pool = prepare_asset_pool() do slot, model
            asset_counter[] += 1
            return (asset_slot=slot, asset_model=String(model.name), asset_token=asset_counter[])
        end
        prepared_assets_model = prepare_model(dummy_prepared_parallel, 0.55f0, 16; name=:asset_model, context=ctx, assets=asset_pool, pool_size=2)
        out_assets, s_assets = prop_run(prepared_assets_model; PASSVALUE=1, slot=1)
        @test size(out_assets) == (16, 16)
        @test s_assets == 1.0f0
        @test asset_counter[] == 1
        out_assets_repeat, _ = prop_run(prepared_assets_model; PASSVALUE=1, slot=1)
        @test size(out_assets_repeat) == (16, 16)
        @test asset_counter[] == 1
        Proper.ensure_prepared_batch_contexts!(prepared_assets_model.batch, 3)
        stack_assets_model, samplings_assets_model = prop_run_multi(prepared_assets_model; PASSVALUE=1:3)
        @test size(stack_assets_model) == (16, 16, 3)
        @test samplings_assets_model == fill(1.0f0, 3)
        @test asset_counter[] == 3
        @test length(unique(vec(stack_assets_model[1, 1, :]))) == 3
        reset_prepared_model!(prepared_assets_model)
        @test all(isnothing, prepared_assets_model.assets.cache)
        _ = prop_run(prepared_assets_model; PASSVALUE=1, slot=1)
        @test asset_counter[] == 4

        prepared_parallel = prepare_prescription(dummy_prepared_parallel, 0.55f0, 16; context=ctx)
        stack_parallel, samplings_parallel = prop_run_multi(prepared_parallel; PASSVALUE=1:3)
        @test size(stack_parallel) == (16, 16, 3)
        @test eltype(stack_parallel) == Float32
        @test length(unique(vec(stack_parallel[1, 1, :]))) == 3
        @test samplings_parallel == fill(1.0f0, 3)

        prepared_parallel_batch = prepare_prescription_batch(prepared_parallel; pool_size=2)
        stack_parallel_batch, samplings_parallel_batch = prop_run_multi(prepared_parallel_batch; PASSVALUE=1:3)
        @test size(stack_parallel_batch) == (16, 16, 3)
        @test eltype(stack_parallel_batch) == Float32
        @test length(unique(vec(stack_parallel_batch[1, 1, :]))) == 3
        @test samplings_parallel_batch == fill(1.0f0, 3)

        prepared_parallel_model = prepare_model(prepared_parallel; name=:dummy_parallel, pool_size=2)
        stack_parallel_model, samplings_parallel_model = prop_run_multi(prepared_parallel_model; PASSVALUE=1:3)
        @test size(stack_parallel_model) == (16, 16, 3)
        @test eltype(stack_parallel_model) == Float32
        @test length(unique(vec(stack_parallel_model[1, 1, :]))) == 3
        @test samplings_parallel_model == fill(1.0f0, 3)

        wf_begin = prop_begin(1.0, 500f-9, 16; context=ctx)
        @test eltype(wf_begin.field) == ComplexF32
        @test wf_begin.workspace === ctx.workspace

        wf_wavefront = prop_wavefront(16, 500f-9, 1.0f0; sampling_m=1f-3, workspace=ctx.workspace)
        @test eltype(wf_wavefront.field) == ComplexF32
        @test wf_wavefront.workspace === ctx.workspace
    end

    @testset "prop_dm application parity" begin
        Random.seed!(1234)
        wf_apply = prop_begin(1.0, 550e-9, 64)
        wf_manual = prop_begin(1.0, 550e-9, 64)
        dm = randn(6, 6) .* 1e-9

        dmap_apply = prop_dm(wf_apply, dm, 2.5, 2.5, 1.0e-3)
        dmap_manual = prop_dm(wf_manual, dm, 2.5, 2.5, 1.0e-3; NO_APPLY=true)

        @test dmap_apply ≈ dmap_manual
        prop_add_phase(wf_manual, 2 .* transpose(dmap_manual))
        @test wf_apply.field ≈ wf_manual.field
    end

    @testset "prop_zernikes application parity" begin
        wf_apply = prop_begin(1.0, 550e-9, 64)
        wf_manual = prop_begin(1.0, 550e-9, 64)

        zmap_apply = prop_zernikes(wf_apply, [4, 7], [20e-9, -15e-9])
        zmap_manual = prop_zernikes(wf_manual, [4, 7], [20e-9, -15e-9]; no_apply=true)

        @test zmap_apply ≈ zmap_manual
        prop_add_phase(wf_manual, zmap_manual)
        @test wf_apply.field ≈ wf_manual.field
    end
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
    @test prop_set_antialiasing(3) == 3
    @test_throws ArgumentError prop_set_antialiasing(4)
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
    pix_out = similar(pix)
    @test prop_pixellate!(pix_out, rand(32, 32), 2) === pix_out

    m1 = prop_polygon(wf, 6, 0.2)
    m2 = prop_irregular_polygon(wf, [-0.1, 0.1, 0.1, -0.1], [-0.1, -0.1, 0.1, 0.1])
    m3 = prop_rounded_rectangle(wf, 0.05, 0.2, 0.3)
    m4 = similar(m3)
    @test prop_rounded_rectangle!(m4, wf, 0.05, 0.2, 0.3) === m4
    @test size(m1) == size(wf.field)
    @test size(m2) == size(wf.field)
    @test size(m3) == size(wf.field)

    @test prop_fftw() == true
    @test prop_ffti() == false
    defs = prop_dftidefs()
    @test defs[:DFTI_FORWARD_SCALE] == 4
    @test defs[:DFTI_BACKWARD_SCALE] == 5
    @test libcconv(rand(8, 8), 4.2, 3.7) isa Real
    @test libcconvthread(rand(8, 8), 4.2, 3.7) isa Real
    @test size(prop_szoom(rand(8, 8), 2.0), 1) == 3
    @test size(prop_szoom(rand(8, 8), 2.0, 16), 1) == 16
    szoom_out = zeros(16, 16)
    @test prop_szoom!(szoom_out, rand(8, 8), 2.0) === szoom_out
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

    # Mutating PSD map path should match wrapper output with same RNG seed.
    dmap_out = zeros(Float64, size(wf.field)...)
    dmap_mut = prop_psd_errormap!(dmap_out, wf, 1e-18, 10.0, 3.0; no_apply=true, rng=MersenneTwister(1))
    @test dmap_mut === dmap_out
    @test dmap_mut == dmap
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
    # Upstream cubic_conv uses 0-based coordinates and edge clamping.
    @test prop_cubic_conv(a, 2.0, 2.0) ≈ a[3, 3] atol=1e-12
    out = prop_cubic_conv(a, [1.0, 2.0, 3.0], [1.0, 2.0]; grid=true)
    @test size(out) == (2, 3)
    @test size(prop_szoom(a, 2.0, 8), 1) == 8
end
