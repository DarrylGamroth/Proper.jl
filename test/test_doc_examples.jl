using Test
using Proper

@testset "Doc examples" begin
    @testset "prop_run example" begin
        function demo_prescription(λm, n)
            wf = prop_begin(1.0, λm, n)
            prop_circular_aperture(wf, 0.5)
            prop_define_entrance(wf)
            return prop_end(wf)
        end

        psf, sampling = prop_run(demo_prescription, 0.55, 32)
        @test size(psf) == (32, 32)
        @test sampling > 0
    end

    @testset "PASSVALUE compatibility shim example" begin
        function pass_demo(λm, n; radius=0.5)
            wf = prop_begin(1.0, λm, n)
            prop_circular_aperture(wf, radius)
            prop_define_entrance(wf)
            return prop_end(wf)
        end

        psf, sampling = prop_run(pass_demo, 0.55, 16; PASSVALUE=Dict(:radius => 0.25))
        @test size(psf) == (16, 16)
        @test sampling > 0
    end

    @testset "prepare_prescription example" begin
        demo_prescription(λm, n) = prop_end(prop_begin(1.0, λm, n))

        prepared = prepare_prescription(demo_prescription, 0.55, 16)
        psf, sampling = prop_run(prepared)
        @test size(psf) == (16, 16)
        @test sampling > 0
    end

    @testset "prepare_prescription explicit precision example" begin
        demo_prescription(λm, n) = prop_end(prop_begin(1.0f0, λm, n))

        prepared = prepare_prescription(demo_prescription, 0.55f0, 16; precision=Float32)
        psf, sampling = prop_run(prepared)
        @test size(psf) == (16, 16)
        @test eltype(psf) == Float32
        @test sampling > 0
    end

    @testset "prepare_prescription_batch example" begin
        function batch_demo(λm, n, pass)
            wf = prop_begin(1.0 + pass, λm, n)
            return prop_end(wf)
        end

        prepared = prepare_prescription(batch_demo, 0.55, 16)
        batch = prepare_prescription_batch(prepared; pool_size=2)
        stack, samplings = prop_run_multi(batch; PASSVALUE=[0.0, 0.5])
        @test size(stack) == (16, 16, 2)
        @test length(samplings) == 2
    end

    @testset "prepare_model example" begin
        function model_demo(λm, n, pass; gain=1.0)
            psf = fill(Float64(gain + pass), n, n)
            return psf, 1.0e-3
        end

        model = prepare_model(:model_demo, model_demo, 0.55, 8; pool_size=2, assets=(gain=2.0,))
        psf, sampling = prop_run(model; slot=1, PASSVALUE=3.0)
        @test psf[1, 1] == 5.0
        @test sampling == 1.0e-3
    end

    @testset "prepare_hot_call example" begin
        function hot_demo(λm, n; payload)
            psf = fill(Float64(payload.gain), n, n)
            return psf, 1.0e-3
        end

        payload = (gain=2.5,)
        model = prepare_model(hot_demo, 0.55, 8; pool_size=1)
        hot = prepare_hot_call(model; payload=payload)
        psf, sampling = prop_run_hot(hot)
        @test psf[1, 1] == 2.5
        @test sampling == 1.0e-3
    end

    @testset "vector of prepared runs example" begin
        sweep_demo(λm, n) = prop_end(prop_begin(1.0f0, λm, n))

        runs = [
            prepare_prescription(sweep_demo, 0.50f0, 8; precision=Float32),
            prepare_prescription(sweep_demo, 0.55f0, 8; precision=Float32),
            prepare_prescription(sweep_demo, 0.60f0, 8; precision=Float32),
        ]

        stack, samplings = prop_run_multi(runs)
        @test size(stack) == (8, 8, 3)
        @test eltype(stack) == Float32
        @test length(samplings) == 3
    end
end
