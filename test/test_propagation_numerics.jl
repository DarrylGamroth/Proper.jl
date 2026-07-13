using LinearAlgebra

function _direct_dft2(a::AbstractMatrix{<:Complex}, direction::Int)
    direction in (-1, 1) || throw(ArgumentError("direction must be -1 or +1"))
    ny, nx = size(a)
    T = float(real(eltype(a)))
    out = Matrix{Complex{T}}(undef, ny, nx)
    two_pi = T(2pi)
    @inbounds for qx in 0:(nx - 1), qy in 0:(ny - 1)
        acc = zero(Complex{T})
        for x in 0:(nx - 1), y in 0:(ny - 1)
            turns = T(qy * y) / T(ny) + T(qx * x) / T(nx)
            acc += a[y + 1, x + 1] * cis(T(direction) * two_pi * turns)
        end
        out[qy + 1, qx + 1] = acc
    end
    return out
end

@inline _oracle_fft_index(k::Int, n::Int) = k <= fld(n - 1, 2) ? k : k - n

function _oracle_spatial_phase(n::Int, dx, coefficient)
    return [
        cis(coefficient * (
            (_oracle_fft_index(i - 1, n) * dx)^2 +
            (_oracle_fft_index(j - 1, n) * dx)^2
        )) for i in 1:n, j in 1:n
    ]
end

function _oracle_frequency_phase(n::Int, dx, coefficient)
    df = inv(n * dx)
    return [
        cis(coefficient * (
            (_oracle_fft_index(i - 1, n) * df)^2 +
            (_oracle_fft_index(j - 1, n) * df)^2
        )) for i in 1:n, j in 1:n
    ]
end

function _unit_contract_prescription(λm, n)
    return fill(λm, n, n), λm
end

@testset "Propagation numerical contract" begin
    @testset "propagation rejects unsupported rectangular grids without mutation" begin
        initial = reshape(ComplexF64.(1:24), 4, 6)
        for (routine, setup!) in (
            (prop_ptp, _ -> nothing),
            (prop_wts, _ -> nothing),
            (prop_stw, wf -> (wf.reference_surface = Proper.SPHERICAL)),
        )
            wf = Proper.WaveFront(copy(initial), 500e-9, 1e-3, 0.0, 1.0)
            setup!(wf)
            field_before = copy(wf.field)
            state_before = (
                wf.z_m,
                wf.sampling_m,
                wf.reference_surface,
                wf.beam_type_old,
                wf.propagator_type,
            )
            err = try
                routine(wf, 0.01)
                nothing
            catch exception
                exception
            end
            @test err isa ArgumentError
            @test occursin("square wavefront field", sprint(showerror, err))
            @test wf.field == field_before
            @test (
                wf.z_m,
                wf.sampling_m,
                wf.reference_surface,
                wf.beam_type_old,
                wf.propagator_type,
            ) == state_before
        end

        wf = Proper.WaveFront(copy(initial), 500e-9, 1e-3, 0.0, 1.0)
        field_before = copy(wf.field)
        state_before = (wf.z_m, wf.beam_type_old, wf.propagator_type)
        @test_throws ArgumentError prop_propagate(wf, 0.01)
        @test wf.field == field_before
        @test (wf.z_m, wf.beam_type_old, wf.propagator_type) == state_before

        wf_zero = Proper.WaveFront(copy(initial), 500e-9, 1e-3, 0.0, 1.0)
        @test_throws ArgumentError prop_ptp(wf_zero, 0.0)
    end

    @testset "FFT propagation matches an independent direct DFT" begin
        n = 5
        λ = 632.8e-9
        dx = 1.7e-3
        dz = 0.037
        initial = reshape(
            ComplexF64.(1:(n * n)) .+ im .* reverse(ComplexF64.(1:(n * n))),
            n,
            n,
        ) ./ 17
        initial_power = sum(abs2, initial)

        wf_ptp = Proper.WaveFront(copy(initial), λ, dx, 0.0, 1.0)
        ptp_spectrum = _direct_dft2(initial, -1) ./ n
        ptp_spectrum .*= _oracle_frequency_phase(n, dx, -pi * λ * dz)
        ptp_expected = _direct_dft2(ptp_spectrum, 1) ./ n
        prop_ptp(wf_ptp, dz)
        @test isapprox(wf_ptp.field, ptp_expected; atol=2e-12, rtol=2e-12)
        @test isapprox(sum(abs2, wf_ptp.field), initial_power; atol=2e-11, rtol=2e-13)

        wf_roundtrip = Proper.WaveFront(copy(initial), λ, dx, 0.0, 1.0)
        prop_ptp(wf_roundtrip, dz)
        prop_ptp(wf_roundtrip, -dz)
        @test isapprox(wf_roundtrip.field, initial; atol=3e-12, rtol=3e-12)
        @test iszero(wf_roundtrip.z_m)

        for d in (dz, -dz)
            wf_wts = Proper.WaveFront(copy(initial), λ, dx, 0.0, 1.0)
            qfield = initial .* _oracle_spatial_phase(n, dx, pi / (λ * d))
            direction = d >= 0 ? -1 : 1
            wts_expected = _direct_dft2(qfield, direction) ./ n
            prop_wts(wf_wts, d)
            @test isapprox(wf_wts.field, wts_expected; atol=2e-12, rtol=2e-12)
            @test isapprox(sum(abs2, wf_wts.field), initial_power; atol=2e-11, rtol=2e-13)
            @test wf_wts.sampling_m == λ * abs(d) / (n * dx)

            wf_stw = Proper.WaveFront(copy(initial), λ, dx, 0.0, 1.0)
            wf_stw.reference_surface = Proper.SPHERICAL
            stw_dx = λ * abs(d) / (n * dx)
            transformed = _direct_dft2(initial, direction) ./ n
            stw_expected = transformed .* _oracle_spatial_phase(n, stw_dx, pi / (λ * d))
            prop_stw(wf_stw, d)
            @test isapprox(wf_stw.field, stw_expected; atol=2e-12, rtol=2e-12)
            @test isapprox(sum(abs2, wf_stw.field), initial_power; atol=2e-11, rtol=2e-13)
            @test wf_stw.sampling_m == stw_dx
        end
    end

    @testset "coordinate and unit conversions are explicit" begin
        @test collect(Proper.coordinate_axis(4, 0.25)) == [-0.5, -0.25, 0.0, 0.25]
        @test collect(Proper.coordinate_axis(5, 0.25)) == [-0.5, -0.25, 0.0, 0.25, 0.5]
        @test collect(Proper.spatial_frequency_axis(4, 0.5)) == [-1.0, -0.5, 0.0, 0.5]
        @test Proper.fft_order_rho2_map(4, 4, 0.5) == [
            0.0 0.25 1.0 0.25
            0.25 0.5 1.25 0.5
            1.0 1.25 2.0 1.25
            0.25 0.5 1.25 0.5
        ]

        wf = prop_begin(2.0, 500e-9, 8)
        wf.current_fratio = 5.0
        wf.sampling_m = 0.01
        @test prop_get_sampling(wf) === 0.01
        @test prop_get_sampling_radians(wf) === 0.001
        @test prop_get_sampling_arcsec(wf) == 0.001 * Proper._RAD_TO_ARCSEC

        wavelength_image, sampling = prop_run(_unit_contract_prescription, 0.55, 3)
        @test all(==(0.55e-6), wavelength_image)
        @test sampling == 0.55e-6

        mktempdir() do dir
            map_value = 0.125
            map_path = joinpath(dir, "unit_map.fits")
            prop_writemap(fill(map_value, 8, 8), map_path; SAMPLING=wf.sampling_m)

            wf_nm = prop_begin(1.0, 500e-9, 8)
            wf_um = prop_begin(1.0, 500e-9, 8)
            prop_errormap(wf_nm, map_path; SAMPLING=wf_nm.sampling_m, NM=true)
            prop_errormap(wf_um, map_path; SAMPLING=wf_um.sampling_m, MICRONS=true)
            @test all(isapprox(cis(2pi * map_value * 1e-9 / wf_nm.wavelength_m); atol=2e-14), wf_nm.field)
            @test all(isapprox(cis(2pi * map_value * 1e-6 / wf_um.wavelength_m); atol=2e-14), wf_um.field)
        end
    end
end
