function _check_wavefront_geometry(geometry, wf)
    for field in fieldnames(typeof(geometry))
        field === :gridsize && continue
        @test isequal(getproperty(geometry, field), getproperty(wf, field))
    end
    @test prop_get_gridsize(geometry) == prop_get_gridsize(wf)
    @test isequal(prop_get_distancetofocus(geometry), prop_get_distancetofocus(wf))
    @test isequal(prop_get_sampling(geometry), prop_get_sampling(wf))
end

@testset "Scalar wavefront geometry" begin
    for T in (Float32, Float64)
        n = 18
        diameter = T(0.002)
        wavelength = T(6.328e-7)
        fraction = T(7 / n)
        new_geometry() = WaveFrontGeometry(diameter, wavelength, n; beam_diam_fraction=fraction)
        new_wavefront() = prop_begin(diameter, wavelength, n; beam_diam_fraction=fraction)

        for direction in (-one(T), one(T))
            wf, geometry = new_wavefront(), new_geometry()
            rayleigh = direction * wf.z_rayleigh_m
            _check_wavefront_geometry(geometry, wf)
            for (distance, transition) in (
                    (rayleigh / 4, Proper.INSIDE_TO_INSIDE),
                    (2rayleigh, Proper.INSIDE_TO_OUTSIDE),
                    (rayleigh / 2, Proper.OUTSIDE_TO_OUTSIDE),
                    (-11rayleigh / 4, Proper.OUTSIDE_TO_INSIDE))
                @test prop_propagate(geometry, distance) === geometry
                prop_propagate(wf, distance)
                @test geometry.propagator_type === transition
                _check_wavefront_geometry(geometry, wf)
            end
        end

        wf, geometry = new_wavefront(), new_geometry()
        for (operation, value) in ((prop_lens, T(0.45)), (prop_propagate, T(0.18)),
                (prop_lens, T(-0.3)), (prop_propagate, T(-0.06)), (prop_lens, T(0.12)))
            @test operation(geometry, value) === geometry
            operation(wf, value)
            _check_wavefront_geometry(geometry, wf)
        end
        prop_propagate(geometry, prop_get_distancetofocus(geometry) + T(0.015); TO_PLANE=true)
        prop_propagate(wf, prop_get_distancetofocus(wf) + T(0.015); TO_PLANE=true)
        _check_wavefront_geometry(geometry, wf)
        @test_throws ArgumentError prop_lens(geometry, zero(T))

        wf, geometry = new_wavefront(), new_geometry()
        for (operation, value) in ((prop_stw, T(0.12)), (prop_wts, -zero(T)),
                (prop_ptp, T(1e-13)), (prop_stw, zero(T)),
                (prop_wts, T(0.12)), (prop_stw, T(-0.12)))
            z_before = geometry.z_m
            operation(geometry, value)
            operation(wf, value)
            _check_wavefront_geometry(geometry, wf)
            if operation === prop_ptp
                @test geometry.z_m === z_before
                @test geometry.reference_surface === Proper.SPHERICAL
                @test_throws ArgumentError prop_ptp(geometry, T(1e-11))
            end
        end
        prop_begin!(wf, diameter, wavelength; beam_diam_fraction=fraction)
        _check_wavefront_geometry(new_geometry(), wf)

        wf, geometry = new_wavefront(), new_geometry()
        for aperture in (prop_circular_aperture, prop_circular_obscuration)
            @test aperture(geometry, diameter / 2) === geometry
            aperture(wf, diameter / 2)
            _check_wavefront_geometry(geometry, wf)
            @test_throws InexactError aperture(wf, floatmax(T))
            @test_throws InexactError aperture(geometry, floatmax(T))
        end
        @test all(isbitstype, fieldtypes(typeof(geometry)))
    end
    @test_throws ArgumentError WaveFrontGeometry(1.0, 1e-6, 0)
end

@testset "Scalar geometry phase validation" begin
    for (diameter, wavelength, short_distance) in (
            (1f-24, 1f-6, 1f-36), (1e-159, 1e-6, 1e-308)),
            make_wavefront in (prop_begin, WaveFrontGeometry)
        T = typeof(wavelength)
        # An inside-to-outside transition applies a spatial phase whose
        # coefficient overflows despite finite positive input parameters.
        wf = make_wavefront(diameter, wavelength, 8)
        @test_throws DomainError prop_propagate(wf, short_distance)

        wf = make_wavefront(diameter, wavelength, 8)
        @test_throws DomainError prop_ptp(wf, one(T))

        wf = make_wavefront(diameter / T(10), wavelength, 8)
        prop_wts(wf, zero(T))
        @test_throws DomainError prop_stw(wf, one(T))

        wf = make_wavefront(T(0.002), wavelength, 8)
        wf.sampling_m = floatmax(T)
        @test_throws DomainError prop_lens(wf, one(T))

        # A one-pixel grid contains only the axis origin. An infinite phase
        # coefficient produces NaN there and the ordinary phase fill accepts it.
        wf = make_wavefront(diameter, wavelength, 1)
        @test prop_propagate(wf, short_distance) === wf
    end
end
