using FFTW

@testset "Public helper and branch coverage" begin
    @testset "small public wavefront helpers" begin
        wf = prop_begin(1.0, 500e-9, 8)
        wf.field .= reshape(ComplexF64.(1:64), 8, 8) .* cis(0.25)

        @test prop_get_wavefront(wf) === wf.field
        @test prop_get_amplitude(wf) ≈ abs.(wf.field)
        @test prop_get_phase(wf) ≈ angle.(wf.field)
        @test prop_get_refradius(wf) == wf.z_m - wf.z_w0_m
        @test prop_get_distancetofocus(wf) == wf.z_w0_m - wf.z_m
        @test prop_get_nyquistsampling(wf) == wf.current_fratio * wf.wavelength_m / 2
        @test prop_get_nyquistsampling(wf, 600e-9) == wf.current_fratio * 600e-9 / 2
        @test prop_get_sampling_arcsec(wf) > 0

        other = fill(ComplexF64(2, 3), 8, 8)
        before = copy(wf.field)
        @test prop_add_wavefront(wf, other) === wf
        @test wf.field == before .+ other
        @test_throws ArgumentError prop_add_wavefront(wf, ones(4, 4))
    end

    @testset "scalar map operations and simple utilities" begin
        wf_mul = prop_begin(1.0, 500e-9, 5)
        prop_multiply(wf_mul, 2)
        @test all(wf_mul.field .== ComplexF64(2))

        wf_div = prop_begin(1.0, 500e-9, 5)
        prop_divide(wf_div, 2)
        @test all(wf_div.field .== ComplexF64(0.5))

        wf_phase = prop_begin(1.0, 500e-9, 5)
        prop_add_phase(wf_phase, 1e-9)
        @test all(wf_phase.field .≈ cis((2pi / wf_phase.wavelength_m) * 1e-9))

        mktempdir() do dir
            state_path = joinpath(dir, "state")
            @test prop_init_savestate(state_path) == state_path
            @test isdir(state_path)

            wisdom_path = joinpath(dir, "wisdom.fftw")
            @test prop_fftw_wisdom(wisdom_path) == wisdom_path
            @test isfile(wisdom_path)
            @test prop_load_fftw_wisdom(wisdom_path)
            @test !prop_load_fftw_wisdom(joinpath(dir, "missing.fftw"))
        end

        grid_wisdom = Proper._wisdom_path(2, 1)
        rm(grid_wisdom; force=true)
        try
            @test !prop_load_fftw_wisdom(3, 1)
            @test prop_fftw_wisdom(2; nthreads=1) == grid_wisdom
            @test isfile(grid_wisdom)
            @test prop_load_fftw_wisdom(2, 1)
        finally
            rm(grid_wisdom; force=true)
        end

        mktempdir() do dir
            outpath = joinpath(dir, "zernikes.txt")
            open(outpath, "w") do io
                redirect_stdout(io) do
                    @test Proper.prop_print_zernikes(3) === nothing
                end
            end
            @test occursin("j=1", read(outpath, String))
        end
    end

    @testset "keyword and grid helpers" begin
        normalized = Proper.normalize_kwargs(pairs((; FOO=1, bar=2)))
        @test normalized.foo == 1
        @test normalized.bar == 2
        @test_throws ArgumentError Proper.normalize_kwargs(pairs(NamedTuple{(:FOO, :foo)}((1, 2))))
        @test Proper.compat_bool(true)
        @test Proper.compat_bool(1)
        @test !Proper.compat_bool(0)
        @test_throws ArgumentError Proper.compat_bool("true")
        @test Proper.kw_resolve(nothing, 2, 3) == 2
        @test Proper.kw_resolve(nothing, nothing, 3) == 3
        @test Proper.kw_resolve_bool(nothing, 1, false)
        @test Proper.kw_resolve_float(nothing, 2, 0.0) == 2.0
        @test Proper.kw_resolve_float(nothing, nothing, nothing) === nothing
        @test Proper.kw_resolve_string(nothing, :abc, nothing) == "abc"
        @test Proper.kw_resolve_symbol(nothing, "METH", :linear) == :meth
        @test Proper.kw_lookup(pairs((; METH="cubic")), :METH, "linear") == "cubic"
        @test Proper.kw_lookup(pairs((; meth="cubic")), :METH, "linear") == "cubic"
        @test Proper.kw_lookup_bool(pairs((; DARK=1)), :DARK, false)
        @test Proper.kw_lookup_float(pairs((; X=2)), :X, 0.0) == 2.0
        @test Proper.kw_lookup_float(pairs((;)), :x, nothing) === nothing
        @test Proper.kw_lookup_string(pairs((; NAME=:demo)), :NAME, nothing) == "demo"
        @test Proper.kw_lookup_present(pairs((; X=1)), :X)
        @test Proper.kw_lookup_present(pairs((; x=1)), :X)

        @test collect(Proper.coordinate_axis(4, 0.5)) == [-1.0, -0.5, 0.0, 0.5]
        r = Proper.radius_map(3, 3, 1.0)
        @test size(r) == (3, 3)
        @test r[2, 2] == 0
        @test collect(Proper.spatial_frequency_axis(4, 0.5)) == [-1.0, -0.5, 0.0, 0.5]
        rho2 = Proper.fft_order_rho2_map(4, 4, 0.5)
        rsqr = Proper.fft_order_rsqr_map(4, 4, 0.5)
        @test size(rho2) == (4, 4)
        @test size(rsqr) == (4, 4)
        @test rho2[1, 1] == 0
        @test minimum(rsqr) == 0
    end

    @testset "geometry and zernike branch coverage" begin
        wf = prop_begin(1.0, 500e-9, 32)
        offaxis = prop_begin(1.0, 500e-9, 32)
        dark = prop_begin(1.0, 500e-9, 32)
        prop_circular_aperture(wf, 0.25; NORM=true)
        prop_circular_aperture(offaxis, 0.25, 0.05, -0.05; NORM=true)
        prop_circular_aperture(dark, 0.25; NORM=true, DARK=true)
        @test any(!iszero, abs.(wf.field))
        @test any(!iszero, abs.(offaxis.field))
        @test any(iszero, abs.(dark.field))
        @test_throws ArgumentError Proper.circle_geometry(Float64, prop_begin(1.0, 500e-9, 16), 0.0, 0.0, 0.0, Proper.CircleOptions(Float64, pairs((;))))

        @test Proper.CircleOptions(pairs((; NORM=true, DARK=true))).norm
        @test Proper.shifted_circle_apply_exec_style(Proper.GeometryKAStyle()) isa Proper.ShiftedCircleKAExecStyle
        centered_geom = Proper.circle_geometry(Float64, wf, 0.25, 0.0, 0.0, Proper.CircleOptions(Float64, pairs((; NORM=true))))
        shifted_geom = Proper.circle_geometry(Float64, wf, 0.25, 0.05, -0.05, Proper.CircleOptions(Float64, pairs((; NORM=true))))
        @test Proper.circle_center_exec_style(centered_geom) isa Proper.CenteredCircleStyle
        @test Proper.circle_center_exec_style(shifted_geom) isa Proper.ShiftedCircleStyle
        @test Proper.circle_mask_value(100.0, 100.0, centered_geom, false, 2) == 0.0
        @test Proper.circle_mask_value(0.0, 0.0, centered_geom, false, 2) == 1.0
        @test Proper.circle_mask_value(centered_geom.rad_pix, 0.0, centered_geom, true, 2) < 1.0

        wf_ka_center = prop_begin(1.0, 500e-9, 32)
        wf_ka_shift = prop_begin(1.0, 500e-9, 32)
        wf_ka_invert = prop_begin(1.0, 500e-9, 32)
        @test Proper._apply_shifted_circle!(Proper.ShiftedCircleKAExecStyle(), wf_ka_center, 0.25, 0.0, 0.0, Proper.CircleOptions(Float64, pairs((; NORM=true))), false) === wf_ka_center
        @test Proper._apply_shifted_circle!(Proper.ShiftedCircleKAExecStyle(), wf_ka_shift, 0.25, 0.05, -0.05, Proper.CircleOptions(Float64, pairs((; NORM=true))), false) === wf_ka_shift
        @test Proper._apply_shifted_circle!(Proper.ShiftedCircleKAExecStyle(), wf_ka_invert, 0.25, 0.0, 0.0, Proper.CircleOptions(Float64, pairs((; NORM=true))), true) === wf_ka_invert
        @test isapprox(wf_ka_center.field, wf.field; atol=1e-12, rtol=1e-12)
        @test isapprox(wf_ka_shift.field, offaxis.field; atol=1e-12, rtol=1e-12)
        @test any(iszero, abs.(wf_ka_invert.field))

        wfz = prop_begin(1.0, 500e-9, 16)
        @test Proper._facti(-1) == 0.0
        maps = prop_zernikes(wfz, 4)
        @test size(maps) == (16, 16, 4)
        dmap = prop_zernikes(wfz, [1, 2], [1e-9, 2e-9]; no_apply=true)
        @test size(dmap) == (16, 16)
        @test all(isfinite, dmap)
        amp_wf = prop_begin(1.0, 500e-9, 16)
        amp_map = prop_zernikes(amp_wf, 1, 0.5; amplitude=true, radius=0.5)
        @test size(amp_map) == (16, 16)
        @test all(isfinite, amp_map)
        obscured = prop_zernikes(prop_begin(1.0, 500e-9, 16), 4, 1e-9, 0.2; no_apply=true)
        @test all(isfinite, obscured)
        @test_throws ArgumentError prop_zernikes(wfz, [1, 2], [1e-9])
        @test_throws ArgumentError prop_zernikes(wfz, 23, 1e-9, 0.1)
    end

    @testset "cubic convolution and propagation branch coverage" begin
        a = reshape(collect(1.0:25.0), 5, 5)
        @test prop_cubic_conv(a, 2.0, 2.0) isa Real
        @test prop_cubic_conv(RunContext(typeof(a)), a, 2.0, 2.0) isa Real
        x = [2.0, 3.0, 4.0]
        y = [2.0, 3.0, 4.0]
        pointwise = prop_cubic_conv(a, x, y; grid=false)
        @test size(pointwise) == (3,)
        @test_throws ArgumentError prop_cubic_conv(a, [2.0, 3.0], [2.0]; grid=false)
        grid = prop_cubic_conv(a, x, y; grid=true)
        @test size(grid) == (3, 3)
        out = similar(grid)
        @test prop_cubic_conv_grid!(out, RunContext(typeof(a)), a, x, y) === out
        @test out ≈ grid
        @test_throws ArgumentError prop_cubic_conv_grid!(zeros(2, 2), a, x, y)
        xgrid = [2.0 3.0; 4.0 2.5]
        ygrid = [2.0 2.5; 3.0 4.0]
        coord = prop_cubic_conv(a, xgrid, ygrid)
        @test size(coord) == size(xgrid)
        @test_throws ArgumentError prop_cubic_conv(a, xgrid, ygrid[1:1, :])

        wf_ptp = prop_begin(1.0, 500e-9, 16)
        wf_ptp.reference_surface = Proper.SPHERICAL
        @test_throws ArgumentError prop_ptp(wf_ptp, 0.01)

        wf_ptp_zero = prop_begin(1.0, 500e-9, 16)
        @test prop_ptp(wf_ptp_zero, 0.0, RunContext(wf_ptp_zero), Proper.FFTWorkspace(Float64)) === wf_ptp_zero

        wf_ptp_generic = prop_begin(1.0, 500e-9, 16)
        ctx_generic = RunContext(wf_ptp_generic)
        @test Proper._prop_ptp_fft!(
            Proper.PTPGenericExecStyle(),
            wf_ptp_generic,
            ctx_generic,
            Proper.fft_workspace(ctx_generic),
            16,
            16,
            16,
            wf_ptp_generic.sampling_m,
            -pi * wf_ptp_generic.wavelength_m * 0.01,
        ) === wf_ptp_generic

        wf_wts = prop_begin(1.0, 500e-9, 16)
        @test prop_wts(wf_wts, -0.01) === wf_wts
        @test wf_wts.reference_surface === Proper.SPHERICAL

        wf_wts_generic_pos = prop_begin(1.0, 500e-9, 16)
        wf_wts_generic_neg = prop_begin(1.0, 500e-9, 16)
        @test Proper._prop_wts_fft_step!(Proper.GenericFFTStyle(), wf_wts_generic_pos, 0.01, 16, RunContext(wf_wts_generic_pos), Proper.FFTWorkspace(Float64)) === wf_wts_generic_pos
        @test Proper._prop_wts_fft_step!(Proper.GenericFFTStyle(), wf_wts_generic_neg, -0.01, 16, RunContext(wf_wts_generic_neg), Proper.FFTWorkspace(Float64)) === wf_wts_generic_neg

        wf_stw = prop_begin(1.0, 500e-9, 16)
        wf_stw.reference_surface = Proper.SPHERICAL
        wf_stw.z_w0_m = 0.02
        @test prop_stw(wf_stw, -0.01) === wf_stw
        @test wf_stw.reference_surface === Proper.PLANAR

        wf_stw_generic_pos = prop_begin(1.0, 500e-9, 16)
        wf_stw_generic_neg = prop_begin(1.0, 500e-9, 16)
        @test Proper._prop_stw_fft_step!(Proper.GenericFFTStyle(), wf_stw_generic_pos, 0.01, 16, RunContext(wf_stw_generic_pos), Proper.FFTWorkspace(Float64)) === wf_stw_generic_pos
        @test Proper._prop_stw_fft_step!(Proper.GenericFFTStyle(), wf_stw_generic_neg, -0.01, 16, RunContext(wf_stw_generic_neg), Proper.FFTWorkspace(Float64)) === wf_stw_generic_neg

        wf_fallback = prop_begin(1.0, 500e-9, 16)
        wf_fallback.reference_surface = Proper.PLANAR
        @test prop_stw(wf_fallback, 0.01) === wf_fallback
    end
end
