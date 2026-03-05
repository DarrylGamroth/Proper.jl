using Test
using Random

@testset "R3 mutating kernels and workspace reuse" begin
    @testset "Mutating interpolation kernels" begin
        rng = MersenneTwister(42)
        a = rand(rng, Float64, 16, 16)
        x = collect(range(0.0, 15.0; length=9))
        y = collect(range(0.0, 15.0; length=7))

        out = similar(a, length(y), length(x))
        Proper.prop_cubic_conv_grid!(out, a, x, y)
        ref = prop_cubic_conv(a, x, y; grid=true)
        @test isapprox(out, ref; atol=0, rtol=0)

        # Warm-up + allocation gate for grid interpolation kernel.
        Proper.prop_cubic_conv_grid!(out, a, x, y)
        @test (@allocated Proper.prop_cubic_conv_grid!(out, a, x, y)) == 0
    end

    @testset "Context workspace reuse" begin
        rng = MersenneTwister(123)
        wf = prop_begin(1.0, 500e-9, 32)
        ctx = RunContext(wf)
        ctx2 = RunContext(wf)
        @test Proper.interp_workspace(ctx) === Proper.interp_workspace(ctx2)
        @test Proper.fft_workspace(ctx) === Proper.fft_workspace(ctx2)
        dmap = rand(rng, Float32, 32, 32)
        opts = Proper.ResampleMapOptions(wf, wf.sampling_m, 16.0, 16.0)
        out = zeros(Float64, size(wf.field)...)

        Proper.prop_resamplemap!(out, wf, dmap, opts, ctx)
        iws = Proper.interp_workspace(ctx)
        xptr = pointer(iws.xcoords)
        yptr = pointer(iws.ycoords)

        Proper.prop_resamplemap!(out, wf, dmap, opts, ctx)
        @test pointer(iws.xcoords) == xptr
        @test pointer(iws.ycoords) == yptr
        @test length(iws.xcoords) == size(wf.field, 2)
        @test length(iws.ycoords) == size(wf.field, 1)

        # Warm-up + allocation gate for context-routed resample path.
        Proper.prop_resamplemap!(out, wf, dmap, opts, ctx)
        @test (@allocated Proper.prop_resamplemap!(out, wf, dmap, opts, ctx)) == 0

        wfptp = prop_begin(1.0, 500e-9, 32)
        Proper.prop_ptp(wfptp, 0.01, ctx)
        fws = Proper.fft_workspace(ctx)
        rhoptr = pointer(fws.rho2)
        Proper.prop_ptp(wfptp, 0.01, ctx)
        @test pointer(fws.rho2) == rhoptr
        @test fws.forward_plan !== nothing
        @test fws.backward_plan !== nothing

        # Wavefront-owned mask buffer is reused by aperture wrappers.
        wfmask = prop_begin(1.0, 500e-9, 32)
        Proper.prop_circular_aperture(wfmask, 0.2)
        mws = wfmask.workspace.mask
        mptr = pointer(mws.mask)
        Proper.prop_circular_aperture(wfmask, 0.2)
        @test pointer(mws.mask) == mptr
        Proper.prop_elliptical_aperture(wfmask, 0.2, 0.15)
        @test pointer(mws.mask) == mptr
        Proper.prop_elliptical_obscuration(wfmask, 0.2, 0.15)
        @test pointer(mws.mask) == mptr
        Proper.prop_rectangular_aperture(wfmask, 0.25, 0.18)
        @test pointer(mws.mask) == mptr
        Proper.prop_rectangular_obscuration(wfmask, 0.25, 0.18)
        @test pointer(mws.mask) == mptr
        Proper.prop_circular_aperture(wfmask, 0.2) # warmup
        @test (@allocated Proper.prop_circular_aperture(wfmask, 0.2)) < 50_000
        Proper.prop_elliptical_aperture(wfmask, 0.2, 0.15) # warmup
        @test (@allocated Proper.prop_elliptical_aperture(wfmask, 0.2, 0.15)) < 70_000
        Proper.prop_rectangular_aperture(wfmask, 0.25, 0.18) # warmup
        @test (@allocated Proper.prop_rectangular_aperture(wfmask, 0.25, 0.18)) < 70_000
    end

    @testset "Mutating wrapper parity" begin
        rng = MersenneTwister(7)
        img = rand(rng, Float32, 16, 16)
        ctx = RunContext(typeof(img))

        rout = similar(img)
        Proper.prop_rotate!(rout, img, 12.0, ctx)
        rref = prop_rotate(img, 12.0, ctx)
        @test isapprox(rout, rref; atol=0, rtol=0)

        mout = similar(img, 16, 16)
        Proper.prop_magnify!(mout, img, 1.2, ctx; QUICK=true)
        mref = prop_magnify(img, 1.2, 16, ctx; QUICK=true)
        @test isapprox(mout, mref; atol=1e-6, rtol=1e-6)
    end

    @testset "Mutating prop_end reuse" begin
        wf = prop_begin(1.0, 500e-9, 32)
        prop_circular_aperture(wf, 0.2)
        prop_propagate(wf, 0.01)

        iref, sref = prop_end(wf)
        ibuf = similar(iref)
        @test (@inferred Proper.prop_end!(ibuf, wf)) === ibuf
        iout, s = prop_end(wf, ibuf)
        @test iout === ibuf
        @test s == sref
        @test isapprox(iout, iref; atol=0, rtol=0)
        Proper.prop_end!(ibuf, wf) # warmup
        @test (@allocated Proper.prop_end!(ibuf, wf)) == 0

        cref, _ = prop_end(wf; noabs=true)
        cbuf = similar(cref)
        @test (@inferred Proper.prop_end!(cbuf, wf; noabs=true)) === cbuf
        @test isapprox(Proper.prop_end!(cbuf, wf; noabs=true), cref; atol=0, rtol=0)
        Proper.prop_end!(cbuf, wf; noabs=true) # warmup
        @test (@allocated Proper.prop_end!(cbuf, wf; noabs=true)) == 0

        xref, _ = prop_end(wf; extract=8)
        xbuf = similar(xref)
        @test (@inferred Proper.prop_end!(xbuf, wf; extract=8)) === xbuf
        @test isapprox(Proper.prop_end!(xbuf, wf; extract=8), xref; atol=0, rtol=0)
        Proper.prop_end!(xbuf, wf; extract=8) # warmup
        @test (@allocated Proper.prop_end!(xbuf, wf; extract=8)) == 0
    end

    @testset "KA pilot parity on large CPU arrays" begin
        n = 512
        @test !Proper.ka_mask_enabled(Matrix{ComplexF64}, n, n)
        @test Proper.ka_end_enabled(Matrix{ComplexF64}, n, n)

        wfmask = prop_begin(1.0, 500e-9, n)
        RT = real(eltype(wfmask.field))
        mask = zeros(RT, n, n)
        Proper.prop_ellipse!(mask, wfmask, 0.35, 0.35, 0.0, 0.0; NORM=true)
        f_ka = copy(wfmask.field)
        f_loop = copy(wfmask.field)
        Proper.ka_apply_shifted_mask!(f_ka, mask)
        Proper._apply_shifted_mask_loop!(f_loop, mask)
        @test isapprox(f_ka, f_loop; atol=0, rtol=0)

        wfend = prop_begin(1.0, 500e-9, n)
        prop_circular_aperture(wfend, 0.2)
        prop_propagate(wfend, 0.01)
        out_ka = Matrix{Float64}(undef, n, n)
        out_loop = similar(out_ka)
        Proper.prop_end!(out_ka, wfend)
        Proper._copy_shifted_intensity_loop!(out_loop, wfend.field, 1, 1)
        @test isapprox(out_ka, out_loop; atol=0, rtol=0)
    end

    @testset "Mutating geometry parity" begin
        wf = prop_begin(1.0, 500e-9, 64)
        RT = real(eltype(wf.field))
        n = size(wf.field, 1)

        ellipse_out = zeros(RT, n, n)
        Proper.prop_ellipse!(ellipse_out, wf, 0.35, 0.25, 0.05, -0.03; ROTATION=13.0, DARK=true, NORM=true)
        ellipse_ref = prop_ellipse(wf, 0.35, 0.25, 0.05, -0.03; ROTATION=13.0, DARK=true, NORM=true)
        @test isapprox(ellipse_out, ellipse_ref; atol=0, rtol=0)

        rect_out = zeros(RT, n, n)
        Proper.prop_rectangle!(rect_out, wf, 0.4, 0.2, 0.03, -0.05; ROTATION=22.0, NORM=true)
        rect_ref = prop_rectangle(wf, 0.4, 0.2, 0.03, -0.05; ROTATION=22.0, NORM=true)
        @test isapprox(rect_out, rect_ref; atol=0, rtol=0)

        poly_out = zeros(RT, n, n)
        Proper.prop_polygon!(poly_out, wf, 6, 0.33, 0.02, -0.04; ROTATION=7.0, NORM=true)
        poly_ref = prop_polygon(wf, 6, 0.33, 0.02, -0.04; ROTATION=7.0, NORM=true)
        @test isapprox(poly_out, poly_ref; atol=0, rtol=0)

        xverts = RT[-0.20, 0.12, 0.28, -0.08]
        yverts = RT[-0.18, -0.22, 0.19, 0.25]
        ipoly_out = zeros(RT, n, n)
        Proper.prop_irregular_polygon!(ipoly_out, wf, xverts, yverts; NORM=true)
        ipoly_ref = prop_irregular_polygon(wf, xverts, yverts; NORM=true)
        @test isapprox(ipoly_out, ipoly_ref; atol=0, rtol=0)
    end
end
