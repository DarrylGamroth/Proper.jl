using Test
using Random

function _matlab_round(x::Real)
    return x > 0 ? floor(x + 0.5) : -floor(-x + 0.5)
end

function _matlab_szoom_ref(arri::AbstractMatrix, magn::Real; nox::Integer=0, noy::Integer=0)
    fsp = 6
    ks = 1 + 2 * fsp
    niy, nix = size(arri)
    nox = nox > 0 ? nox : floor(Int, nix * magn) - ks
    noy = noy > 0 ? noy : nox
    iicx = (nix ÷ 2) + 1
    iocx = (nox ÷ 2) + 1
    iicy = (niy ÷ 2) + 1
    iocy = (noy ÷ 2) + 1

    function build_table(nout)
        tbl = Matrix{Float64}(undef, nout, ks)
        vk = collect(0:(ks - 1)) .- (ks ÷ 2)
        vp = (collect(0:(nout - 1)) .- (nout ÷ 2)) ./ magn
        vp .-= map(_matlab_round, vp)
        for row in 1:nout
            for col in 1:ks
                at = (vk[col] - vp[row]) * pi
                if abs(vk[col] - vp[row]) <= fsp
                    tbl[row, col] = at == 0 ? 1.0 : (sin(at) / at) * (sin(at / fsp) / (at / fsp))
                else
                    tbl[row, col] = 0.0
                end
            end
        end
        return tbl
    end

    tblx = build_table(nox)
    tbly = noy == nox ? tblx : build_table(noy)
    out = zeros(Float64, noy, nox)
    for ioy in 1:noy
        iiy = Int(_matlab_round((ioy - iocy) / magn) + iicy)
        iiy1 = iiy - fsp
        iiy2 = iiy + fsp
        if iiy1 < 1 || iiy2 > niy
            continue
        end
        strp = arri[iiy1:iiy2, :] .* tbly[ioy, :]
        for iox in 1:nox
            iix = Int(_matlab_round((iox - iocx) / magn) + iicx)
            iix1 = iix - fsp
            iix2 = iix + fsp
            if iix1 < 1 || iix2 > nix
                continue
            end
            out[ioy, iox] = sum(sum(strp[:, iix1:iix2], dims=1)[:] .* tblx[iox, :])
        end
    end
    return out
end

function _matlab_pixellate_ref(aimi::AbstractMatrix, sami::Real, samo::Real, npo::Integer=0)
    npiy, npix = size(aimi)
    icix = (npix ÷ 2) + 1
    iciy = (npiy ÷ 2) + 1
    magn = sami / samo
    vx = Proper.FFTW.ifftshift(collect(1:npix) .- icix) ./ (icix - 1) ./ 2.0 ./ magn
    vy = Proper.FFTW.ifftshift(collect(1:npiy) .- iciy) ./ (iciy - 1) ./ 2.0 ./ magn
    pmtf = [sinc(y) * sinc(x) for y in vy, x in vx]
    amtf = pmtf .* Proper.fft(Proper.FFTW.ifftshift(aimi))
    cimc = Proper.FFTW.fftshift(abs.(Proper.ifft(amtf)) ./ (magn * magn))
    out_n = npo > 0 ? npo : floor(Int, magn * npix)
    return prop_magnify(cimc, magn, out_n)
end

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

    @testset "prop_resamplemap axis construction matches MATLAB formula" begin
        wf = prop_begin(1.0, 500e-9, 5)
        dmap = reshape(collect(1.0:25.0), 5, 5)
        pix = wf.sampling_m

        # MATLAB:
        # o1x = ([1:nx] - fix(nx/2) - 1) * bm.dx / pxsc + cpx - xshift/pxsc
        # o1y = ([1:ny] - fix(ny/2) - 1) * bm.dx / pxsc + cpy - yshift/pxsc
        xshift = pix
        yshift = -pix
        xcoords = ((collect(1.0:5.0) .- floor(5 / 2) .- 1) .* (wf.sampling_m / pix)) .+ 2.0 .- xshift / pix
        ycoords = ((collect(1.0:5.0) .- floor(5 / 2) .- 1) .* (wf.sampling_m / pix)) .+ 2.0 .- yshift / pix

        ref = Proper.prop_cubic_conv_grid!(Matrix{Float64}(undef, 5, 5), dmap, xcoords, ycoords)
        got = prop_resamplemap(wf, dmap, pix, 2.0, 2.0, xshift, yshift)
        @test got == ref
    end

    @testset "prop_magnify QUICK axis construction matches MATLAB formula" begin
        a = reshape(collect(1.0:16.0), 4, 4)
        mag = 2.0
        nox = 8
        noy = 8
        cx = floor(size(a, 2) / 2) + 1
        cy = floor(size(a, 1) / 2) + 1
        sx = floor(nox / 2) + 1
        sy = floor(noy / 2) + 1
        xcoords = ((collect(1.0:nox) .- sx) ./ mag) .+ cx
        ycoords = ((collect(1.0:noy) .- sy) ./ mag) .+ cy

        ref = Proper.prop_cubic_conv_grid!(Matrix{Float64}(undef, noy, nox), a, xcoords, ycoords)
        got = prop_magnify(a, mag, nox; QUICK=true)
        @test got == ref
    end

    @testset "prop_magnify default output sizing uses fix semantics" begin
        a = reshape(collect(1.0:64.0), 8, 8)
        @test size(prop_magnify(a, 1.6), 1) == 12
        @test size(prop_magnify(a, 1.6; QUICK=true), 1) == 12
    end

    @testset "prop_szoom matches MATLAB rounding semantics" begin
        a = reshape(collect(1.0:64.0), 8, 8)
        ref = _matlab_szoom_ref(a, 0.95; nox=16, noy=16)
        got = prop_szoom(a, 0.95, 16)
        @test isapprox(got, ref; atol=1e-12, rtol=1e-12)
    end

    @testset "prop_pixellate sampling overload matches MATLAB formula" begin
        a = reshape(collect(1.0:64.0), 8, 8)
        ref = _matlab_pixellate_ref(a, 0.5, 1.0, 4)
        got = prop_pixellate(a, 0.5, 1.0, 4)
        @test isapprox(got, ref; atol=1e-12, rtol=1e-12)
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

    @testset "prop_psd_errormap amplitude map is normalized even with NO_APPLY" begin
        wf = prop_begin(1.0, 500e-9, 32)
        field0 = copy(wf.field)
        amp_target = 0.9
        dmap = prop_psd_errormap(wf, 1e-18, 10.0, 3.0; AMPLITUDE=amp_target, no_apply=true, rng=MersenneTwister(9))

        @test maximum(dmap) ≈ amp_target atol=1e-12 rtol=0
        @test wf.field == field0
    end

    @testset "DM correction behavior matches upstream trend" begin
        exdir = joinpath(@__DIR__, "..", "examples")
        mod = load_example_module(joinpath(exdir, "run_coronagraph_dm.jl"))
        run_coronagraph_dm = getfield(mod, :run_coronagraph_dm)

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
