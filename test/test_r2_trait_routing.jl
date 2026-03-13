using Test

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
        dmap = rand(Float32, 16, 16)
        opts = Proper.ResampleMapOptions(wf, wf.sampling_m, 8.0, 8.0)
        ctx = RunContext(typeof(dmap))

        out_ctx = similar(dmap, Float64, size(wf.field)...)
        out_def = similar(dmap, Float64, size(wf.field)...)
        Proper.prop_resamplemap!(out_ctx, wf, dmap, opts, ctx)
        Proper.prop_resamplemap!(out_def, wf, dmap, opts)
        @test isapprox(out_ctx, out_def; atol=0, rtol=0)

        img = rand(Float32, 16, 16)
        r_ctx = prop_rotate(img, 12.0, ctx)
        r_def = prop_rotate(img, 12.0)
        @test isapprox(r_ctx, r_def; atol=0, rtol=0)

        m_ctx = prop_magnify(img, 1.2, 16, ctx; QUICK=true)
        m_def = prop_magnify(img, 1.2, 16; QUICK=true)
        @test isapprox(m_ctx, m_def; atol=1e-6, rtol=1e-6)
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

        img = rand(Float32, n, n)
        pix_loop = prop_pixellate(img, 2)
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
    end

    @testset "Context-routed propagation kernels" begin
        wf1 = prop_begin(1.0, 500e-9, 32)
        wf2 = prop_begin(1.0, 500e-9, 32)
        ctx = RunContext(typeof(wf2.field))

        prop_propagate(wf1, 0.05)
        prop_propagate(wf2, 0.05, ctx)

        @test wf1.z_m == wf2.z_m
        @test isapprox(wf1.field, wf2.field; atol=0, rtol=0)
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
            p = prop_pixellate(a, 2)
            @test size(m) == (16, 16)
            @test size(r) == (16, 16)
            @test size(s) == (16, 16)
            @test size(p) == (8, 8)

            wf = Proper.WaveFront(CUDA.fill(ComplexF32(1), 16, 16), 500f-9, 1f-3, 0f0, 1f0)
            prop_qphase(wf, 0.25f0, ctx)
            wf.reference_surface = Proper.PLANAR
            prop_ptp(wf, 0.01f0, ctx)
            prop_circular_aperture(wf, 2.5f-4)
            rect = prop_rectangle(wf, 5f-4, 4f-4)
            round = prop_rounded_rectangle(wf, 2f-4, 5f-4, 4f-4)
            out, sampling = prop_end(wf)
            @test size(rect) == (16, 16)
            @test size(round) == (16, 16)
            @test size(out) == (16, 16)
            @test sampling == wf.sampling_m
        else
            @test true
        end
    end
end
