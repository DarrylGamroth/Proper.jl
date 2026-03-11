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
            a = CUDA.rand(Float32, 16, 16)
            ctx = RunContext(typeof(a))
            @test ctx.backend isa Proper.BackendStyle
            @test ctx.interp isa Proper.InterpStyle

            m = prop_magnify(a, 1.1, 16, ctx; QUICK=true)
            r = prop_rotate(a, 5.0, ctx)
            @test size(m) == (16, 16)
            @test size(r) == (16, 16)
        else
            @test true
        end
    end
end
