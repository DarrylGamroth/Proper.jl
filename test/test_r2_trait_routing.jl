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
        @test isapprox(m_ctx, m_def; atol=0, rtol=0)
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
