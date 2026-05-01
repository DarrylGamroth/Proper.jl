using Test
using Random

function _warmed_gpu_qphase_alloc(wf, z, ctx, sync!)
    prop_qphase(wf, z, ctx)
    sync!()
    prop_qphase(wf, z, ctx)
    sync!()
    return @allocated begin
        prop_qphase(wf, z, ctx)
        sync!()
    end
end

function _warmed_gpu_ptp_alloc(wf, dz, ctx, sync!)
    wf.reference_surface = Proper.PLANAR
    prop_ptp(wf, dz, ctx)
    sync!()
    wf.reference_surface = Proper.PLANAR
    prop_ptp(wf, dz, ctx)
    sync!()
    return @allocated begin
        wf.reference_surface = Proper.PLANAR
        prop_ptp(wf, dz, ctx)
        sync!()
    end
end

function _warmed_gpu_wts_alloc(wf, dz, ctx, sync!)
    wf.reference_surface = Proper.PLANAR
    prop_wts(wf, dz, ctx)
    sync!()
    wf.reference_surface = Proper.PLANAR
    prop_wts(wf, dz, ctx)
    sync!()
    return @allocated begin
        wf.reference_surface = Proper.PLANAR
        prop_wts(wf, dz, ctx)
        sync!()
    end
end

function _warmed_gpu_stw_alloc(wf, dz, ctx, sync!)
    wf.reference_surface = Proper.SPHERICAL
    prop_stw(wf, dz, ctx)
    sync!()
    wf.reference_surface = Proper.SPHERICAL
    prop_stw(wf, dz, ctx)
    sync!()
    return @allocated begin
        wf.reference_surface = Proper.SPHERICAL
        prop_stw(wf, dz, ctx)
        sync!()
    end
end

function _warmed_gpu_end_real_alloc(out, wf, sync!)
    wf.reference_surface = Proper.PLANAR
    Proper.prop_end!(out, wf)
    sync!()
    wf.reference_surface = Proper.PLANAR
    Proper.prop_end!(out, wf)
    sync!()
    return @allocated begin
        wf.reference_surface = Proper.PLANAR
        Proper.prop_end!(out, wf)
        sync!()
    end
end

function _warmed_gpu_end_complex_alloc(out, wf, sync!)
    wf.reference_surface = Proper.PLANAR
    Proper.prop_end!(out, wf; noabs=true)
    sync!()
    wf.reference_surface = Proper.PLANAR
    Proper.prop_end!(out, wf; noabs=true)
    sync!()
    return @allocated begin
        wf.reference_surface = Proper.PLANAR
        Proper.prop_end!(out, wf; noabs=true)
        sync!()
    end
end

const GPU_WARM_QPHASE_ALLOC_MAX = 12_288
const GPU_WARM_PTP_ALLOC_MAX = 12_288
const GPU_WARM_WTS_ALLOC_MAX = 12_288
const GPU_WARM_STW_ALLOC_MAX = 8_192
const GPU_WARM_END_REAL_ALLOC_MAX = 16_384
const GPU_WARM_END_COMPLEX_ALLOC_MAX = 16_384

function _gpu_map_apply_smoke!(
    wf_gpu,
    gpu_real_matrix,
    expected_backend_array_type,
    sync!,
)
    n = size(wf_gpu.field, 1)
    scale_map_cpu = fill(Float32(1.25), n, n)
    phase_map_cpu = fill(Float32(1f-9), n, n)
    scale_map_gpu = copyto!(similar(gpu_real_matrix, size(scale_map_cpu)...), scale_map_cpu)
    phase_map_gpu = copyto!(similar(gpu_real_matrix, size(phase_map_cpu)...), phase_map_cpu)

    wf_mul_gpu = Proper.WaveFront(similar(wf_gpu.field), wf_gpu.wavelength_m, wf_gpu.sampling_m, wf_gpu.z_m, wf_gpu.beam_diameter_m)
    fill!(wf_mul_gpu.field, ComplexF32(1))
    wf_mul_cpu = Proper.WaveFront(fill(ComplexF32(1), n, n), wf_gpu.wavelength_m, wf_gpu.sampling_m, wf_gpu.z_m, wf_gpu.beam_diameter_m)
    prop_multiply(wf_mul_gpu, scale_map_gpu)
    sync!()
    prop_multiply(wf_mul_cpu, scale_map_cpu)
    @test isapprox(Array(wf_mul_gpu.field), wf_mul_cpu.field; atol=1f-6, rtol=1f-6)

    wf_div_gpu = Proper.WaveFront(similar(wf_gpu.field), wf_gpu.wavelength_m, wf_gpu.sampling_m, wf_gpu.z_m, wf_gpu.beam_diameter_m)
    fill!(wf_div_gpu.field, ComplexF32(2))
    wf_div_cpu = Proper.WaveFront(fill(ComplexF32(2), n, n), wf_gpu.wavelength_m, wf_gpu.sampling_m, wf_gpu.z_m, wf_gpu.beam_diameter_m)
    prop_divide(wf_div_gpu, scale_map_gpu)
    sync!()
    prop_divide(wf_div_cpu, scale_map_cpu)
    @test isapprox(Array(wf_div_gpu.field), wf_div_cpu.field; atol=1f-6, rtol=1f-6)

    wf_phase_gpu = Proper.WaveFront(similar(wf_gpu.field), wf_gpu.wavelength_m, wf_gpu.sampling_m, wf_gpu.z_m, wf_gpu.beam_diameter_m)
    fill!(wf_phase_gpu.field, ComplexF32(1))
    wf_phase_cpu = Proper.WaveFront(fill(ComplexF32(1), n, n), wf_gpu.wavelength_m, wf_gpu.sampling_m, wf_gpu.z_m, wf_gpu.beam_diameter_m)
    prop_add_phase(wf_phase_gpu, phase_map_gpu)
    sync!()
    prop_add_phase(wf_phase_cpu, phase_map_cpu)
    @test isapprox(Array(wf_phase_gpu.field), wf_phase_cpu.field; atol=1f-6, rtol=1f-6)

    mktempdir() do dir
        f = joinpath(dir, "map.fits")
        Proper.prop_fits_write(f, scale_map_cpu; HEADER=Dict("PIXSIZE" => wf_gpu.sampling_m))

        map_read = prop_readmap(wf_gpu, f; SAMPLING=wf_gpu.sampling_m)
        sync!()
        @test map_read isa expected_backend_array_type

        wf_err_gpu = Proper.WaveFront(similar(wf_gpu.field), wf_gpu.wavelength_m, wf_gpu.sampling_m, wf_gpu.z_m, wf_gpu.beam_diameter_m)
        fill!(wf_err_gpu.field, ComplexF32(1))
        wf_err_cpu = Proper.WaveFront(fill(ComplexF32(1), n, n), wf_gpu.wavelength_m, wf_gpu.sampling_m, wf_gpu.z_m, wf_gpu.beam_diameter_m)
        prop_errormap(wf_err_gpu, f; SAMPLING=wf_gpu.sampling_m, WAVEFRONT=true)
        sync!()
        prop_errormap(wf_err_cpu, f; SAMPLING=wf_gpu.sampling_m, WAVEFRONT=true)
        @test isapprox(Array(wf_err_gpu.field), wf_err_cpu.field; atol=1f-5, rtol=1f-5)
    end

    wf_psd_gpu = Proper.WaveFront(similar(wf_gpu.field), wf_gpu.wavelength_m, wf_gpu.sampling_m, wf_gpu.z_m, wf_gpu.beam_diameter_m)
    fill!(wf_psd_gpu.field, ComplexF32(1))
    dmap_gpu = prop_psd_errormap(wf_psd_gpu, 1e-18, 10.0, 3.0; no_apply=true, rng=Random.MersenneTwister(7))
    sync!()
    @test dmap_gpu isa expected_backend_array_type

    wf_psd_apply_gpu = Proper.WaveFront(similar(wf_gpu.field), wf_gpu.wavelength_m, wf_gpu.sampling_m, wf_gpu.z_m, wf_gpu.beam_diameter_m)
    fill!(wf_psd_apply_gpu.field, ComplexF32(1))
    wf_psd_apply_cpu = Proper.WaveFront(fill(ComplexF32(1), n, n), wf_gpu.wavelength_m, wf_gpu.sampling_m, wf_gpu.z_m, wf_gpu.beam_diameter_m)
    prop_psd_errormap(wf_psd_apply_gpu, 1e-18, 10.0, 3.0; MIRROR=true, rng=Random.MersenneTwister(11))
    sync!()
    prop_psd_errormap(wf_psd_apply_cpu, 1e-18, 10.0, 3.0; MIRROR=true, rng=Random.MersenneTwister(11))
    @test isapprox(Array(wf_psd_apply_gpu.field), wf_psd_apply_cpu.field; atol=3f-4, rtol=1f-3)
end

@testset "R2 trait-driven routing" begin
    @testset "Style dispatch for cubic convolution" begin
        struct TestInterpStyle <: Proper.InterpStyle end
        Proper._cubic_sample(::TestInterpStyle, a::AbstractMatrix, y::Real, x::Real) = 123.0

        a = reshape(collect(1.0:16.0), 4, 4)
        out = Proper.prop_cubic_conv(TestInterpStyle(), a, [1.0, 2.0], [1.0, 2.0]; grid=true)
        @test size(out) == (2, 2)
        @test all(out .== 123.0)
    end

    @testset "Context-routed interpolation kernels" begin
        wf = prop_begin(1.0, 500e-9, 16)
        dmap = rand(TEST_RNG, Float32, 16, 16)
        opts = Proper.ResampleMapOptions(wf, wf.sampling_m, 8.0, 8.0)
        ctx = RunContext(typeof(dmap))

        out_ctx = similar(dmap, Float64, size(wf.field)...)
        out_def = similar(dmap, Float64, size(wf.field)...)
        Proper.prop_resamplemap!(out_ctx, wf, dmap, opts, ctx)
        Proper.prop_resamplemap!(out_def, wf, dmap, opts)
        @test isapprox(out_ctx, out_def; atol=0, rtol=0)

        img = rand(TEST_RNG, Float32, 16, 16)
        r_ctx = prop_rotate(img, 12.0, ctx)
        r_def = prop_rotate(img, 12.0)
        @test isapprox(r_ctx, r_def; atol=0, rtol=0)

        m_ctx = prop_magnify(img, 1.2, 16, ctx; QUICK=true)
        m_def = prop_magnify(img, 1.2, 16; QUICK=true)
        @test isapprox(m_ctx, m_def; atol=1e-6, rtol=1e-6)
    end

    @testset "Explicit backend boundary for hot output paths" begin
        wf = prop_begin(1.0, 500e-9, 16)
        cpu_out = zeros(Float64, 16, 16)
        @test prop_end!(cpu_out, wf) === cpu_out
    end

    @testset "KA interpolation pilot parity on large CPU arrays" begin
        n = 256
        @test !Proper.ka_cubic_grid_enabled(Matrix{Float32}, n, n)
        @test !Proper.ka_rotate_enabled(Matrix{Float32}, n, n)

        a = reshape(collect(Float32, 1:(n * n)), n, n)
        x = collect(Float32, 1:n)
        y = collect(Float32, 1:n)
        out_loop = similar(a)
        out_ka = similar(a)
        Proper._prop_cubic_conv_grid_loop!(out_loop, Proper.CubicInterpStyle(), a, x, y)
        Proper.ka_cubic_conv_grid!(out_ka, a, x, y)
        @test isapprox(out_ka, out_loop; atol=0, rtol=0)

        opts_cubic = Proper.RotateOptions(a, pairs((; CUBIC=true)))
        rc_loop = similar(a)
        rc_ka = similar(a)
        c = cos(-deg2rad(9.0))
        s = sin(-deg2rad(9.0))
        Proper._prop_rotate_cubic!(Proper.CubicInterpStyle(), rc_loop, a, c, s, opts_cubic)
        Proper.ka_rotate_cubic!(rc_ka, a, c, s, opts_cubic.cx, opts_cubic.cy, opts_cubic.sx, opts_cubic.sy)
        @test isapprox(rc_ka, rc_loop; atol=0, rtol=0)

        opts_linear = Proper.RotateOptions(a, pairs((; METH="linear")))
        rl_loop = similar(a)
        rl_ka = similar(a)
        Proper._prop_rotate_linear!(rl_loop, a, c, s, opts_linear)
        Proper.ka_rotate_linear!(rl_ka, a, c, s, opts_linear.cx, opts_linear.cy, opts_linear.sx, opts_linear.sy)
        @test isapprox(rl_ka, rl_loop; atol=0, rtol=0)
    end

    @testset "KA geometry and sampling parity on large CPU arrays" begin
        n = 256
        @test !Proper.ka_geometry_enabled(Matrix{Float32}, n, n)
        @test !Proper.ka_sampling_enabled(Matrix{Float32}, n, n)

        wf = prop_begin(1.0, 500e-9, n)
        RT = real(eltype(wf.field))

        rect_loop = zeros(RT, n, n)
        rect_ka = similar(rect_loop)
        Proper.prop_rectangle!(rect_loop, wf, 0.4, 0.2, 0.03, -0.05; ROTATION=22.0, NORM=true)
        dx = RT(prop_get_sampling(wf))
        beamrad = RT(prop_get_beamradius(wf))
        pr = beamrad / dx
        θ = RT(deg2rad(22.0))
        cθ = cos(θ)
        sθ = sin(θ)
        xcp = RT(n ÷ 2) + RT(0.03) * pr
        ycp = RT(n ÷ 2) - RT(0.05) * pr
        xrp = RT(0.4) * pr / RT(2)
        yrp = RT(0.2) * pr / RT(2)
        xp = (-xrp, -xrp, xrp, xrp)
        yp = (-yrp, yrp, yrp, -yrp)
        xbox = ntuple(i -> xp[i] * cθ - yp[i] * sθ + xcp, 4)
        ybox = ntuple(i -> xp[i] * sθ + yp[i] * cθ + ycp, 4)
        minx = max(0, floor(Int, minimum(xbox) - one(RT)))
        maxx = min(n - 1, ceil(Int, maximum(xbox) + one(RT)))
        miny = max(0, floor(Int, minimum(ybox) - one(RT)))
        maxy = min(n - 1, ceil(Int, maximum(ybox) + one(RT)))
        Proper.ka_rectangle_mask!(rect_ka, xcp, ycp, xrp, yrp, cθ, sθ, minx, maxx, miny, maxy; nsub=Proper.antialias_subsampling())
        @test isapprox(rect_ka, rect_loop; atol=1e-12, rtol=1e-12)

        ellipse_loop = zeros(RT, n, n)
        ellipse_ka = similar(ellipse_loop)
        Proper.prop_ellipse!(ellipse_loop, wf, 0.35, 0.25, 0.05, -0.03; ROTATION=13.0, DARK=true, NORM=true)
        beamrad_pix = RT(prop_get_beamradius(wf)) / dx
        xcenter_pix = RT(n ÷ 2) + RT(0.05) * beamrad_pix
        ycenter_pix = RT(n ÷ 2) - RT(0.03) * beamrad_pix
        xrad_pix = RT(0.35) * beamrad_pix
        yrad_pix = RT(0.25) * beamrad_pix
        t = RT(deg2rad(13.0))
        sint = sin(t)
        cost = cos(t)
        delx = inv(xrad_pix)
        dely = inv(yrad_pix)
        drx = delx * cost - dely * sint
        dry = delx * sint + dely * cost
        dr = max(abs(drx), abs(dry))
        Proper.ka_ellipse_mask!(ellipse_ka, xcenter_pix, ycenter_pix, xrad_pix, yrad_pix, sint, cost, one(RT) + dr, one(RT) - dr, RT(1 + 1e-10); dark=true, nsub=Proper.antialias_subsampling())
        @test isapprox(ellipse_ka, ellipse_loop; atol=1e-12, rtol=1e-12)

        xverts = RT[-0.20, 0.12, 0.28, -0.08]
        yverts = RT[-0.18, -0.22, 0.19, 0.25]
        ipoly_loop = zeros(RT, n, n)
        ipoly_ka = similar(ipoly_loop)
        Proper.prop_irregular_polygon!(ipoly_loop, wf, xverts, yverts; NORM=true)
        xv = copy(xverts) .* beamrad
        yv = copy(yverts) .* beamrad
        Proper.ka_irregular_polygon_mask!(ipoly_ka, xv, yv, n ÷ 2, n ÷ 2, dx; nsub=Proper.antialias_subsampling())
        @test isapprox(ipoly_ka, ipoly_loop; atol=1e-12, rtol=1e-12)

        round_loop = Proper.prop_rounded_rectangle(wf, 0.05, 0.3, 0.2, 0.01, -0.02)
        round_ka = similar(round_loop)
        Proper.ka_rounded_rectangle_mask!(round_ka, dx, RT(0.01), RT(-0.02), RT(0.05), RT(0.3) / RT(2), RT(0.2) / RT(2))
        @test isapprox(round_ka, round_loop; atol=1e-12, rtol=1e-12)

        img = rand(TEST_RNG, Float32, n, n)
        pix_loop = Proper._prop_pixellate_factor(img, 2)
        pix_ka = similar(pix_loop)
        Proper.ka_pixellate!(pix_ka, img, 2)
        @test isapprox(pix_ka, pix_loop; atol=0, rtol=0)

        mag = Float32(1.35)
        n_out = 128
        szoom_loop = prop_szoom(img, mag, n_out)
        szoom_ka = similar(szoom_loop)
        table_loop = Matrix{Float32}(undef, n_out, Proper.SZOOM_K)
        table_ka = similar(table_loop)
        Proper._fill_szoom_table_loop!(table_loop, mag)
        Proper.ka_szoom_table!(table_ka, mag, n_out, Proper.SZOOM_K, Proper.SZOOM_DK)
        @test isapprox(table_ka, table_loop; atol=0, rtol=0)
        Proper.ka_szoom_apply!(szoom_ka, img, table_ka, mag)
        @test isapprox(szoom_ka, szoom_loop; atol=1e-6, rtol=1e-6)

        rect_img = rand(TEST_RNG, Float32, n, n + 16)
        rect_loop = prop_szoom(rect_img, mag; NOX=120, NOY=96)
        rect_ka = similar(rect_loop)
        tablex_loop = Matrix{Float32}(undef, 120, Proper.SZOOM_K)
        tabley_loop = Matrix{Float32}(undef, 96, Proper.SZOOM_K)
        tablex_ka = similar(tablex_loop)
        tabley_ka = similar(tabley_loop)
        Proper._fill_szoom_table_loop!(tablex_loop, mag)
        Proper._fill_szoom_table_loop!(tabley_loop, mag)
        Proper.ka_szoom_table!(tablex_ka, mag, 120, Proper.SZOOM_K, Proper.SZOOM_DK)
        Proper.ka_szoom_table!(tabley_ka, mag, 96, Proper.SZOOM_K, Proper.SZOOM_DK)
        @test isapprox(tablex_ka, tablex_loop; atol=0, rtol=0)
        @test isapprox(tabley_ka, tabley_loop; atol=0, rtol=0)
        Proper.ka_szoom_apply!(rect_ka, rect_img, tablex_ka, tabley_ka, mag)
        @test isapprox(rect_ka, rect_loop; atol=1e-6, rtol=1e-6)
    end

    @testset "Context-routed propagation kernels" begin
        wf1 = prop_begin(1.0, 500e-9, 32)
        wf2 = prop_begin(1.0, 500e-9, 32)
        ctx = RunContext(typeof(wf2.field))

        prop_propagate(wf1, 0.05)
        prop_propagate(wf2, 0.05, ctx)

        @test wf1.z_m == wf2.z_m
        @test isapprox(wf1.field, wf2.field; atol=0, rtol=0)

        lens1 = prop_begin(1.0, 500e-9, 32)
        lens2 = prop_begin(1.0, 500e-9, 32)
        lensctx = RunContext(lens2)
        prop_lens(lens1, 0.35)
        prop_lens(lens2, 0.35, lensctx)
        @test lens1.reference_surface === lens2.reference_surface
        @test lens1.beam_type_old === lens2.beam_type_old
        @test lens1.propagator_type === lens2.propagator_type
        @test isapprox(lens1.field, lens2.field; atol=0, rtol=0)
    end

    @testset "Optional CUDA smoke (no scalar indexing)" begin
        cuda_ready = false
        try
            @eval using CUDA
            cuda_ready = CUDA.functional()
        catch
            cuda_ready = false
        end

        if cuda_ready
            CUDA.allowscalar(false)
            @test Base.get_extension(Proper, :ProperCUDAExt) !== nothing

            a = CUDA.rand(Float32, 16, 16)
            ctx = RunContext(typeof(a))
            @test ctx.backend isa Proper.CUDABackend
            @test ctx.fft isa Proper.CUFFTStyle
            @test ctx.interp isa Proper.CubicInterpStyle
            @test Proper.ka_cubic_grid_enabled(typeof(a), 16, 16)
            @test Proper.ka_rotate_enabled(typeof(a), 16, 16)
            @test Proper.ka_geometry_enabled(typeof(a), 16, 16)
            @test Proper.ka_sampling_enabled(typeof(a), 16, 16)

            m = prop_magnify(a, 1.1, 16, ctx; QUICK=true)
            r = prop_rotate(a, 5.0, ctx)
            s = prop_szoom(a, 1.1, 16)
            p = Proper._prop_pixellate_factor(a, 2)
            wf_resample = Proper.WaveFront(CUDA.fill(ComplexF32(1), 16, 16), 500f-9, 1f-3, 0f0, 1f0)
            ropts = Proper.ResampleMapOptions(wf_resample, wf_resample.sampling_m, 8f0, 8f0)
            res = similar(a)
            Proper.prop_resamplemap!(res, wf_resample, a, ropts, ctx)
            @test_throws ArgumentError Proper.prop_cubic_conv(a, 1.5f0, 1.5f0)
            img_cpu = reshape(Float32.(1:256), 16, 16)
            xcoords_cpu = repeat(reshape(collect(Float32, range(3, 13; length=8)), 1, :), 8, 1)
            ycoords_cpu = repeat(reshape(collect(Float32, range(3, 13; length=8)), :, 1), 1, 8)
            coord_ref = Proper.prop_cubic_conv(img_cpu, xcoords_cpu, ycoords_cpu; grid=false)
            coord_gpu = Proper.prop_cubic_conv(CUDA.CuArray(img_cpu), xcoords_cpu, ycoords_cpu; grid=false)
            @test coord_gpu isa CUDA.CuArray
            @test isapprox(Array(coord_gpu), coord_ref; atol=1f-5, rtol=1f-5)
            @test size(m) == (16, 16)
            @test size(r) == (16, 16)
            @test size(s) == (16, 16)
            @test size(p) == (8, 8)
            @test size(res) == (16, 16)
            @test Proper.interp_workspace(ctx).xcoords isa CUDA.CuArray
            @test Proper.interp_workspace(ctx).ycoords isa CUDA.CuArray

            wf = Proper.WaveFront(CUDA.fill(ComplexF32(1), 16, 16), 500f-9, 1f-3, 0f0, 1f0)
            prop_qphase(wf, 0.25f0, ctx)
            wf.reference_surface = Proper.PLANAR
            prop_ptp(wf, 0.01f0, ctx)
            fws = Proper.fft_workspace(ctx)
            @test fws.scratch isa CUDA.CuArray
            rho2 = Proper.ensure_rho2_map!(fws, 16, 16, 1f-3)
            @test rho2 isa CUDA.CuArray
            @test fws.forward_plan !== nothing
            @test fws.backward_plan !== nothing
            pfft = fws.forward_plan
            pbfft = fws.backward_plan
            wf.reference_surface = Proper.PLANAR
            prop_ptp(wf, 0.01f0, ctx)
            @test fws.forward_plan === pfft
            @test fws.backward_plan === pbfft
            @test Proper.ensure_rho2_map!(fws, 16, 16, 1f-3) === rho2
            prop_circular_aperture(wf, 2.5f-4)
            @test wf.workspace.mask.mask isa CUDA.CuArray
            wf_ref = Proper.WaveFront(fill(ComplexF32(1), 16, 16), 500f-9, 1f-3, 0f0, 1f0)
            prop_qphase(wf_ref, 0.25f0)
            wf_ref.reference_surface = Proper.PLANAR
            prop_ptp(wf_ref, 0.01f0)
            prop_circular_aperture(wf_ref, 2.5f-4)
            @test isapprox(Array(wf.field), wf_ref.field; atol=1f-5, rtol=1f-5)
            rect = prop_rectangle(wf, 5f-4, 4f-4)
            round = prop_rounded_rectangle(wf, 2f-4, 5f-4, 4f-4)
            out, sampling = prop_end(wf)
            @test_throws ArgumentError Proper.prop_end!(zeros(Float32, 16, 16), wf)
            out_ref, sampling_ref = prop_end(wf_ref)
            @test size(rect) == (16, 16)
            @test size(round) == (16, 16)
            @test size(out) == (16, 16)
            @test sampling == wf.sampling_m
            @test sampling == sampling_ref
            @test isapprox(Array(out), out_ref; atol=1f-5, rtol=1f-5)
            wf_alloc = Proper.WaveFront(CUDA.fill(ComplexF32(1), 16, 16), 500f-9, 1f-3, 0f0, 1f0)
            out_alloc = similar(wf_alloc.field, Float32, 16, 16)
            @test _warmed_gpu_qphase_alloc(wf_alloc, 0.25f0, ctx, CUDA.synchronize) <= GPU_WARM_QPHASE_ALLOC_MAX
            @test _warmed_gpu_ptp_alloc(wf_alloc, 0.01f0, ctx, CUDA.synchronize) <= GPU_WARM_PTP_ALLOC_MAX
            @test _warmed_gpu_wts_alloc(wf_alloc, 0.01f0, ctx, CUDA.synchronize) <= GPU_WARM_WTS_ALLOC_MAX
            @test _warmed_gpu_stw_alloc(wf_alloc, 0.01f0, ctx, CUDA.synchronize) <= GPU_WARM_STW_ALLOC_MAX
            @test _warmed_gpu_end_real_alloc(out_alloc, wf_alloc, CUDA.synchronize) <= GPU_WARM_END_REAL_ALLOC_MAX
            @test _warmed_gpu_end_complex_alloc(similar(wf_alloc.field), wf_alloc, CUDA.synchronize) <= GPU_WARM_END_COMPLEX_ALLOC_MAX
            _gpu_map_apply_smoke!(wf_alloc, CUDA.zeros(Float32, 16, 16), CUDA.CuArray, CUDA.synchronize)
        else
            @test true
        end
    end

    @testset "Optional AMDGPU smoke (no scalar indexing)" begin
        amdgpu_ready = false
        try
            @eval using AMDGPU
            amdgpu_ready = AMDGPU.functional() && AMDGPU.functional(:rocfft)
        catch
            amdgpu_ready = false
        end

        if amdgpu_ready
            AMDGPU.allowscalar(false)
            @test Base.get_extension(Proper, :ProperAMDGPUExt) !== nothing

            a = AMDGPU.rand(Float32, 16, 16)
            ctx = RunContext(typeof(a))
            @test ctx.backend isa Proper.AMDGPUBackend
            @test ctx.fft isa Proper.ROCFFTStyle
            @test ctx.interp isa Proper.CubicInterpStyle
            @test Proper.ka_cubic_grid_enabled(typeof(a), 16, 16)
            @test Proper.ka_geometry_enabled(typeof(a), 16, 16)
            @test Proper.ka_sampling_enabled(typeof(a), 16, 16)

            p = Proper._prop_pixellate_factor(a, 2)
            @test_throws ArgumentError Proper.prop_cubic_conv(a, 1.5f0, 1.5f0)
            img_cpu = reshape(Float32.(1:256), 16, 16)
            xcoords_cpu = repeat(reshape(collect(Float32, range(3, 13; length=8)), 1, :), 8, 1)
            ycoords_cpu = repeat(reshape(collect(Float32, range(3, 13; length=8)), :, 1), 1, 8)
            coord_ref = Proper.prop_cubic_conv(img_cpu, xcoords_cpu, ycoords_cpu; grid=false)
            coord_gpu = Proper.prop_cubic_conv(AMDGPU.ROCArray(img_cpu), xcoords_cpu, ycoords_cpu; grid=false)
            @test coord_gpu isa AMDGPU.ROCArray
            @test isapprox(Array(coord_gpu), coord_ref; atol=3f-4, rtol=1f-3)
            @test size(p) == (8, 8)

            wf = Proper.WaveFront(AMDGPU.fill(ComplexF32(1), 16, 16), 500f-9, 1f-3, 0f0, 1f0)
            prop_qphase(wf, 0.25f0, ctx)
            wf.reference_surface = Proper.PLANAR
            prop_ptp(wf, 0.01f0, ctx)
            fws = Proper.fft_workspace(ctx)
            @test fws.scratch isa AMDGPU.ROCArray
            rho2 = Proper.ensure_rho2_map!(fws, 16, 16, 1f-3)
            @test rho2 isa AMDGPU.ROCArray
            @test fws.forward_plan !== nothing
            @test fws.backward_plan !== nothing
            pfft = fws.forward_plan
            pbfft = fws.backward_plan
            wf.reference_surface = Proper.PLANAR
            prop_ptp(wf, 0.01f0, ctx)
            @test fws.forward_plan === pfft
            @test fws.backward_plan === pbfft
            @test Proper.ensure_rho2_map!(fws, 16, 16, 1f-3) === rho2
            wf_ref = Proper.WaveFront(fill(ComplexF32(1), 16, 16), 500f-9, 1f-3, 0f0, 1f0)
            prop_qphase(wf_ref, 0.25f0)
            wf_ref.reference_surface = Proper.PLANAR
            prop_ptp(wf_ref, 0.01f0)
            @test isapprox(Array(wf.field), wf_ref.field; atol=3f-4, rtol=1f-3)
            prop_circular_aperture(wf, 2.5f-4)
            @test wf.workspace.mask.mask isa AMDGPU.ROCArray
            prop_circular_aperture(wf_ref, 2.5f-4)
            out, sampling = prop_end(wf)
            @test_throws ArgumentError Proper.prop_end!(zeros(Float32, 16, 16), wf)
            out_ref, sampling_ref = prop_end(wf_ref)
            @test size(out) == (16, 16)
            @test sampling == wf.sampling_m
            @test sampling == sampling_ref
            @test Array(out) isa Matrix{Float32}
            wf_alloc = Proper.WaveFront(AMDGPU.fill(ComplexF32(1), 16, 16), 500f-9, 1f-3, 0f0, 1f0)
            out_alloc = similar(wf_alloc.field, Float32, 16, 16)
            @test _warmed_gpu_qphase_alloc(wf_alloc, 0.25f0, ctx, AMDGPU.synchronize) <= GPU_WARM_QPHASE_ALLOC_MAX
            @test _warmed_gpu_ptp_alloc(wf_alloc, 0.01f0, ctx, AMDGPU.synchronize) <= GPU_WARM_PTP_ALLOC_MAX
            @test _warmed_gpu_wts_alloc(wf_alloc, 0.01f0, ctx, AMDGPU.synchronize) <= GPU_WARM_WTS_ALLOC_MAX
            @test _warmed_gpu_stw_alloc(wf_alloc, 0.01f0, ctx, AMDGPU.synchronize) <= GPU_WARM_STW_ALLOC_MAX
            @test _warmed_gpu_end_real_alloc(out_alloc, wf_alloc, AMDGPU.synchronize) <= GPU_WARM_END_REAL_ALLOC_MAX
            @test _warmed_gpu_end_complex_alloc(similar(wf_alloc.field), wf_alloc, AMDGPU.synchronize) <= GPU_WARM_END_COMPLEX_ALLOC_MAX
        else
            @test true
        end
    end
end
