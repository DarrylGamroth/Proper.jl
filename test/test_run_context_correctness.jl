using FFTW
using Random
using Statistics
using Test

function _carrier_strict_prescription(λm, n)
    wf = prop_begin(1.0, λm, n)
    prop_ptp(wf, λm / 4)
    return copy(wf.field), wf.sampling_m
end

function _carrier_pass_prescription(λm, n, quarter_waves)
    wf = prop_begin(1.0, λm, n)
    prop_ptp(wf, quarter_waves * λm / 4)
    return copy(wf.field), wf.sampling_m
end

function _phase_payload_prescription(λm, n; phase_offset)
    value = phase_offset ? one(λm) : zero(λm)
    return fill(complex(value), n, n), λm
end

function _context_psd_prescription(λm, n, _)
    wf = prop_begin(1.0, λm, n)
    dmap = prop_psd_errormap(wf, 1e-18, 10.0, 3.0; no_apply=true)
    return dmap, wf.sampling_m
end

@inline function _expected_carrier_factor(dz, λ)
    return cispi(2 * rem(dz, λ, RoundNearest) / λ)
end

@testset "RunContext numerical correctness" begin
    @testset "FFTW planning policy preserves numerical results" begin
        rng = MersenneTwister(0xfeed)
        λ = 500e-9
        n = 19
        initial = randn(rng, ComplexF64, n, n)

        for (name, propagate!, dz) in (
            (:ptp_positive, prop_ptp, 0.011),
            (:ptp_negative, prop_ptp, -0.013),
            (:wts_positive, prop_wts, 0.017),
            (:wts_negative, prop_wts, -0.019),
        )
            @testset "$name" begin
                estimate_wf = prop_begin(1.0, λ, n)
                measure_wf = prop_begin(1.0, λ, n)
                estimate_wf.field .= initial
                measure_wf.field .= initial
                estimate_ctx = RunContext(
                    typeof(estimate_wf.field),
                    estimate_wf.workspace;
                    fft_planning=Proper.FFTEstimateStyle(),
                )
                measure_ctx = RunContext(
                    typeof(measure_wf.field),
                    measure_wf.workspace;
                    fft_planning=Proper.FFTMeasureStyle(),
                )

                propagate!(estimate_wf, dz, estimate_ctx)
                FFTW.forget_wisdom()
                propagate!(measure_wf, dz, measure_ctx)

                @test sum(abs2, measure_wf.field) > 0
                @test isapprox(measure_wf.field, estimate_wf.field; atol=5e-12, rtol=5e-12)
                @test Proper.fft_workspace(measure_ctx).plan_flags == FFTW.MEASURE
            end
        end

        for dz in (0.023, -0.029)
            estimate_wf = prop_begin(1.0, λ, n)
            measure_wf = prop_begin(1.0, λ, n)
            estimate_wf.field .= initial
            measure_wf.field .= initial
            estimate_wf.reference_surface = Proper.SPHERICAL
            measure_wf.reference_surface = Proper.SPHERICAL
            estimate_ctx = RunContext(
                typeof(estimate_wf.field),
                estimate_wf.workspace;
                fft_planning=Proper.FFTEstimateStyle(),
            )
            measure_ctx = RunContext(
                typeof(measure_wf.field),
                measure_wf.workspace;
                fft_planning=Proper.FFTMeasureStyle(),
            )

            prop_stw(estimate_wf, dz, estimate_ctx)
            FFTW.forget_wisdom()
            prop_stw(measure_wf, dz, measure_ctx)

            @test sum(abs2, measure_wf.field) > 0
            @test isapprox(measure_wf.field, estimate_wf.field; atol=5e-12, rtol=5e-12)
            @test Proper.fft_workspace(measure_ctx).plan_flags == FFTW.MEASURE
        end

        initial32 = randn(rng, ComplexF32, 18, 18)
        estimate32 = prop_begin(1.0f0, 500f-9, 18)
        measure32 = prop_begin(1.0f0, 500f-9, 18)
        estimate32.field .= initial32
        measure32.field .= initial32
        estimate32_ctx = RunContext(
            typeof(estimate32.field),
            estimate32.workspace;
            fft_planning=Proper.FFTEstimateStyle(),
        )
        measure32_ctx = RunContext(
            typeof(measure32.field),
            measure32.workspace;
            fft_planning=Proper.FFTMeasureStyle(),
        )
        prop_ptp(estimate32, 0.01f0, estimate32_ctx)
        FFTW.forget_wisdom()
        prop_ptp(measure32, 0.01f0, measure32_ctx)
        @test sum(abs2, measure32.field) > 0
        @test isapprox(measure32.field, estimate32.field; atol=2f-5, rtol=2f-5)
    end

    @testset "Active context resolution is workspace-safe" begin
        wf = prop_begin(1.0, 500e-9, 16)
        active = RunContext(
            typeof(wf.field),
            wf.workspace;
            rng=MersenneTwister(11),
            verbose=true,
            fft_planning=Proper.FFTMeasureStyle(),
            carrier_phase=TrackCarrierPhase(),
        )
        unrelated_wf = prop_begin(1.0, 500e-9, 16)

        Proper.with_run_context(active) do
            @test Proper.resolve_run_context(wf) === active
            fallback = Proper.resolve_run_context(unrelated_wf)
            @test fallback !== active
            @test fallback.workspace === unrelated_wf.workspace
            @test Proper.carrier_phase_style(fallback) isa EnvelopeOnly
        end

        reference = prop_begin(1.0, 500e-9, 16)
        reference.field .= wf.field
        reference_ctx = RunContext(
            typeof(reference.field),
            reference.workspace;
            fft_planning=Proper.FFTEstimateStyle(),
            carrier_phase=TrackCarrierPhase(),
        )
        Proper.with_run_context(active) do
            prop_ptp(wf, 0.01)
        end
        prop_ptp(reference, 0.01, reference_ctx)
        @test Proper.fft_workspace(active).plan_flags == FFTW.MEASURE
        @test isapprox(wf.field, reference.field; atol=5e-12, rtol=5e-12)

        child = Proper.fresh_context(active; rng=MersenneTwister(12))
        @test Proper.carrier_phase_style(child) isa TrackCarrierPhase
        @test Proper.fft_planning_style(child) isa Proper.FFTMeasureStyle
        @test child.verbose

        mktempdir() do dir
            map_data = reshape(collect(1.0:256.0), 16, 16) .* 1e-9
            map_file = joinpath(dir, "active_context_map.fits")
            prop_writemap(map_data, map_file; SAMPLING=wf.sampling_m)

            read_active_wf = prop_begin(1.0, 500e-9, 16)
            read_explicit_wf = prop_begin(1.0, 500e-9, 16)
            read_active_ctx = RunContext(
                typeof(read_active_wf.field),
                read_active_wf.workspace;
                fft_planning=Proper.FFTMeasureStyle(),
            )
            read_explicit_ctx = RunContext(
                typeof(read_explicit_wf.field),
                read_explicit_wf.workspace;
                fft_planning=Proper.FFTMeasureStyle(),
            )
            active_map = Proper.with_run_context(read_active_ctx) do
                prop_readmap(read_active_wf, map_file; SAMPLING=read_active_wf.sampling_m)
            end
            explicit_map = prop_readmap(
                read_explicit_wf,
                map_file,
                read_explicit_ctx;
                SAMPLING=read_explicit_wf.sampling_m,
            )
            @test active_map == explicit_map

            error_active_wf = prop_begin(1.0, 500e-9, 16)
            error_explicit_wf = prop_begin(1.0, 500e-9, 16)
            error_active_ctx = RunContext(
                typeof(error_active_wf.field),
                error_active_wf.workspace,
            )
            error_explicit_ctx = RunContext(
                typeof(error_explicit_wf.field),
                error_explicit_wf.workspace,
            )
            Proper.with_run_context(error_active_ctx) do
                prop_errormap(error_active_wf, map_file; SAMPLING=error_active_wf.sampling_m)
            end
            prop_errormap(
                error_explicit_wf,
                map_file,
                error_explicit_ctx;
                SAMPLING=error_explicit_wf.sampling_m,
            )
            @test error_active_wf.field == error_explicit_wf.field
        end
    end

    @testset "Carrier phase tracks coherent optical path" begin
        λ = 500e-9
        n = 16
        rng = MersenneTwister(0xc011)
        initial = randn(rng, ComplexF64, n, n)

        for dz in (λ / 4, -λ / 4)
            envelope_wf = prop_begin(1.0, λ, n)
            carrier_wf = prop_begin(1.0, λ, n)
            envelope_wf.field .= initial
            carrier_wf.field .= initial
            envelope_ctx = RunContext(
                typeof(envelope_wf.field),
                envelope_wf.workspace;
                carrier_phase=EnvelopeOnly(),
            )
            carrier_ctx = RunContext(
                typeof(carrier_wf.field),
                carrier_wf.workspace;
                carrier_phase=TrackCarrierPhase(),
            )
            prop_ptp(envelope_wf, dz, envelope_ctx)
            prop_ptp(carrier_wf, dz, carrier_ctx)
            factor = _expected_carrier_factor(dz, λ)
            @test isapprox(carrier_wf.field, envelope_wf.field .* factor; atol=2e-12, rtol=2e-12)
            @test isapprox(abs2.(carrier_wf.field), abs2.(envelope_wf.field); atol=2e-12, rtol=2e-12)
        end

        for dz in (0.01 + λ / 4, -0.01 - λ / 4)
            envelope_wf = prop_begin(1.0, λ, n)
            carrier_wf = prop_begin(1.0, λ, n)
            envelope_wf.field .= initial
            carrier_wf.field .= initial
            envelope_ctx = RunContext(
                typeof(envelope_wf.field),
                envelope_wf.workspace;
                carrier_phase=EnvelopeOnly(),
            )
            carrier_ctx = RunContext(
                typeof(carrier_wf.field),
                carrier_wf.workspace;
                carrier_phase=TrackCarrierPhase(),
            )
            prop_wts(envelope_wf, dz, envelope_ctx)
            prop_wts(carrier_wf, dz, carrier_ctx)
            factor = _expected_carrier_factor(dz, λ)
            @test isapprox(carrier_wf.field, envelope_wf.field .* factor; atol=2e-12, rtol=2e-12)
        end

        for (stw_arg, setup_distance) in ((0.006 + λ / 4, 0.01), (0.0, 0.01 + λ / 4))
            envelope_wf = prop_begin(1.0, λ, n)
            carrier_wf = prop_begin(1.0, λ, n)
            envelope_wf.field .= initial
            carrier_wf.field .= initial
            envelope_ctx = RunContext(typeof(envelope_wf.field), envelope_wf.workspace)
            carrier_setup_ctx = RunContext(typeof(carrier_wf.field), carrier_wf.workspace)
            prop_wts(envelope_wf, setup_distance, envelope_ctx)
            prop_wts(carrier_wf, setup_distance, carrier_setup_ctx)
            actual_distance = iszero(stw_arg) ? carrier_wf.z_w0_m - carrier_wf.z_m : stw_arg
            carrier_ctx = RunContext(
                typeof(carrier_wf.field),
                carrier_wf.workspace;
                carrier_phase=TrackCarrierPhase(),
            )
            prop_stw(envelope_wf, stw_arg, envelope_ctx)
            prop_stw(carrier_wf, stw_arg, carrier_ctx)
            factor = _expected_carrier_factor(actual_distance, λ)
            @test isapprox(carrier_wf.field, envelope_wf.field .* factor; atol=2e-12, rtol=2e-12)
        end

        envelope_wf = prop_begin(1.0, λ, n)
        carrier_wf = prop_begin(1.0, λ, n)
        envelope_wf.field .= initial
        carrier_wf.field .= initial
        envelope_ctx = RunContext(typeof(envelope_wf.field), envelope_wf.workspace)
        carrier_ctx = RunContext(
            typeof(carrier_wf.field),
            carrier_wf.workspace;
            carrier_phase=TrackCarrierPhase(),
        )
        prop_lens(envelope_wf, 0.25, envelope_ctx)
        prop_lens(carrier_wf, 0.25, carrier_ctx)
        total_distance = 0.7 + λ / 4
        prop_propagate(envelope_wf, total_distance, envelope_ctx)
        prop_propagate(carrier_wf, total_distance, carrier_ctx)
        factor = _expected_carrier_factor(total_distance, λ)
        @test isapprox(carrier_wf.field, envelope_wf.field .* factor; atol=5e-12, rtol=5e-12)

        arm1 = prop_begin(1.0, λ, 8)
        arm2 = prop_begin(1.0, λ, 8)
        arm1_ctx = RunContext(
            typeof(arm1.field),
            arm1.workspace;
            carrier_phase=TrackCarrierPhase(),
        )
        arm2_ctx = RunContext(
            typeof(arm2.field),
            arm2.workspace;
            carrier_phase=TrackCarrierPhase(),
        )
        prop_ptp(arm1, λ, arm1_ctx)
        prop_ptp(arm2, 3λ / 2, arm2_ctx)
        @test maximum(abs2, arm1.field .+ arm2.field) < 1e-28

        envelope1 = prop_begin(1.0, λ, 8)
        envelope2 = prop_begin(1.0, λ, 8)
        prop_ptp(envelope1, λ)
        prop_ptp(envelope2, 3λ / 2)
        @test isapprox(mean(abs2, envelope1.field .+ envelope2.field), 4.0; atol=2e-12, rtol=0)
    end

    @testset "Carrier run options are reserved and typed" begin
        raw_upper, _ = prop_run(_carrier_strict_prescription, 0.5, 8; PHASE_OFFSET=true)
        raw_lower, _ = prop_run(_carrier_strict_prescription, 0.5, 8; phase_offset=1)
        raw_equal, _ = prop_run(
            _carrier_strict_prescription,
            0.5,
            8;
            PHASE_OFFSET=true,
            phase_offset=1,
        )
        raw_default, _ = prop_run(_carrier_strict_prescription, 0.5, 8)
        @test all(isapprox.(raw_upper, im; atol=2e-12, rtol=0))
        @test raw_lower == raw_upper
        @test raw_equal == raw_upper
        @test all(isapprox.(raw_default, 1; atol=2e-12, rtol=0))
        @test_throws ArgumentError prop_run(
            _carrier_strict_prescription,
            0.5,
            8;
            PHASE_OFFSET=true,
            phase_offset=false,
        )

        prepared = prepare_prescription(
            _carrier_strict_prescription,
            0.5,
            8;
            PHASE_OFFSET=true,
        )
        @test Proper.carrier_phase_style(prepared.context) isa TrackCarrierPhase
        @test !haskey(prepared.kwargs, :PHASE_OFFSET)
        @test !haskey(prepared.kwargs, :phase_offset)
        prepared_out, _ = prop_run(prepared)
        overridden_out, _ = prop_run(prepared; phase_offset=false)
        @test all(isapprox.(prepared_out, im; atol=2e-12, rtol=0))
        @test all(isapprox.(overridden_out, 1; atol=2e-12, rtol=0))
        @test Proper.carrier_phase_style(prepared.context) isa TrackCarrierPhase
        @test_throws ArgumentError prop_run(prepared; PHASE_OFFSET=true, phase_offset=false)

        default_prepared = prepare_prescription(_carrier_strict_prescription, 0.5, 8)
        prepared_run = prepare_run(default_prepared; phase_offset=true)
        @test Proper.carrier_phase_style(prepared_context(prepared_run)) isa TrackCarrierPhase
        @test !haskey(prepared_run.kwargs, :phase_offset)
        run_out, _ = prop_run(prepared_run)
        @test all(isapprox.(run_out, im; atol=2e-12, rtol=0))
        @test_throws ArgumentError prepare_run(
            default_prepared;
            PHASE_OFFSET=true,
            phase_offset=false,
        )

        outer = RunContext(Matrix{Float64}; carrier_phase=TrackCarrierPhase())
        Proper.with_run_context(outer) do
            nested_disabled, _ = prop_run(
                _carrier_strict_prescription,
                0.5,
                8;
                phase_offset=false,
            )
            @test all(isapprox.(nested_disabled, 1; atol=2e-12, rtol=0))
            @test Proper.active_run_context() === outer
        end

        payload, _ = prop_run(
            _phase_payload_prescription,
            0.5,
            4;
            PASSVALUE=(phase_offset=true,),
        )
        @test all(==(1 + 0im), payload)

        multipass = prepare_prescription(
            _carrier_pass_prescription,
            0.5,
            8;
            phase_offset=true,
        )
        passes = [1, 2, 3]
        stack, _ = prop_run_multi(multipass; PASSVALUE=passes)
        for (i, quarter_waves) in pairs(passes)
            expected = cispi(quarter_waves / 2)
            @test all(isapprox.(selectdim(stack, 3, i), expected; atol=2e-12, rtol=0))
        end
        batch = prepare_prescription_batch(multipass; pool_size=2)
        model = prepare_model(multipass; pool_size=2)
        @test all(ctx -> Proper.carrier_phase_style(ctx) isa TrackCarrierPhase, batch.contexts)
        @test all(ctx -> Proper.carrier_phase_style(ctx) isa TrackCarrierPhase, model.batch.contexts)
        batch_stack, _ = prop_run_multi(batch; PASSVALUE=passes)
        model_stack, _ = prop_run_multi(model; PASSVALUE=passes)
        @test batch_stack == stack
        @test model_stack == stack
    end

    @testset "Prepared PSD maps consume context RNG and planning" begin
        seed = 0x1234
        context_rng = MersenneTwister(seed)
        expected_rng = MersenneTwister(seed)
        wf = prop_begin(1.0, 500e-9, 16)
        ctx = RunContext(typeof(wf.field), wf.workspace; rng=context_rng)
        actual = Proper.with_run_context(ctx) do
            prop_psd_errormap(wf, 1e-18, 10.0, 3.0; no_apply=true)
        end
        expected_wf = prop_begin(1.0, 500e-9, 16)
        expected = prop_psd_errormap(
            expected_wf,
            1e-18,
            10.0,
            3.0;
            no_apply=true,
            rng=expected_rng,
        )
        @test actual == expected
        @test rand(context_rng, UInt64) == rand(expected_rng, UInt64)

        untouched_seed = 0x5678
        override_seed = 0x9abc
        untouched_rng = MersenneTwister(untouched_seed)
        untouched_reference = MersenneTwister(untouched_seed)
        override_rng = MersenneTwister(override_seed)
        override_reference = MersenneTwister(override_seed)
        override_wf = prop_begin(1.0, 500e-9, 16)
        override_ctx = RunContext(
            typeof(override_wf.field),
            override_wf.workspace;
            rng=untouched_rng,
        )
        override_actual = Proper.with_run_context(override_ctx) do
            prop_psd_errormap(
                override_wf,
                1e-18,
                10.0,
                3.0;
                NO_APPLY=true,
                RNG=override_rng,
            )
        end
        override_expected_wf = prop_begin(1.0, 500e-9, 16)
        override_expected = prop_psd_errormap(
            override_expected_wf,
            1e-18,
            10.0,
            3.0;
            no_apply=true,
            rng=override_reference,
        )
        @test override_actual == override_expected
        @test rand(untouched_rng, UInt64) == rand(untouched_reference, UInt64)

        estimate_wf = prop_begin(1.0, 500e-9, 17)
        measure_wf = prop_begin(1.0, 500e-9, 17)
        estimate_ctx = RunContext(
            typeof(estimate_wf.field),
            estimate_wf.workspace;
            rng=MersenneTwister(77),
            fft_planning=Proper.FFTEstimateStyle(),
        )
        measure_ctx = RunContext(
            typeof(measure_wf.field),
            measure_wf.workspace;
            rng=MersenneTwister(77),
            fft_planning=Proper.FFTMeasureStyle(),
        )
        estimate_map = Proper.with_run_context(estimate_ctx) do
            prop_psd_errormap(estimate_wf, 1e-18, 10.0, 3.0; no_apply=true)
        end
        FFTW.forget_wisdom()
        measure_map = Proper.with_run_context(measure_ctx) do
            prop_psd_errormap(measure_wf, 1e-18, 10.0, 3.0; no_apply=true)
        end
        @test isapprox(measure_map, estimate_map; atol=2e-20, rtol=2e-12)
        @test Proper.fft_workspace(measure_ctx).plan_flags == FFTW.MEASURE

        root1 = RunContext(Matrix{Float64}; rng=MersenneTwister(0xbadc0de))
        root2 = RunContext(Matrix{Float64}; rng=MersenneTwister(0xbadc0de))
        prepared1 = prepare_prescription(
            _context_psd_prescription,
            0.5,
            16;
            context=root1,
        )
        prepared2 = prepare_prescription(
            _context_psd_prescription,
            0.5,
            16;
            context=root2,
        )
        batch1 = prepare_prescription_batch(prepared1; pool_size=4)
        batch2 = prepare_prescription_batch(prepared2; pool_size=4)
        passes = collect(1:4)
        threaded_stack, threaded_sampling = prop_run_multi(batch1; PASSVALUE=passes)
        serial_stack = similar(threaded_stack)
        serial_sampling = similar(threaded_sampling)
        for i in eachindex(passes)
            out, sampling = prop_run(batch2; PASSVALUE=passes[i], slot=i)
            copyto!(selectdim(serial_stack, 3, i), out)
            serial_sampling[i] = sampling
        end
        @test threaded_stack == serial_stack
        @test threaded_sampling == serial_sampling
        @test all(i -> !isequal(selectdim(threaded_stack, 3, 1), selectdim(threaded_stack, 3, i)), 2:4)
    end
end
