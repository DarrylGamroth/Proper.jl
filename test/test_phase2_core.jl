using Test

@testset "Phase 2 core kernels" begin
    wf = prop_begin(1.0, 500e-9, 32)
    orig = copy(wf.field)

    prop_lens(wf, 10.0)
    @test wf.field != orig

    wf2 = prop_begin(1.0, 500e-9, 32)
    prop_circular_aperture(wf2, 0.2)
    @test any(iszero, abs.(wf2.field))

    wf3 = prop_begin(1.0, 500e-9, 32)
    prop_propagate(wf3, 0.25)
    @test isapprox(prop_get_z(wf3), 0.25; atol=1e-12)

    a = reshape(collect(1.0:16.0), 4, 4)
    @test size(prop_magnify(a, 2.0)) == (8, 8)
    @test size(prop_rotate(a, 15.0)) == size(a)
end

@testset "Phase 2 run entrypoints" begin
    dummy(λm, n; kwargs...) = prop_begin(1.0, λm, n)

    out, s = prop_run(dummy, 0.55, 16)
    @test size(out) == (16, 16)
    @test s > 0

    passvals = [1, 2, 3]
    stack, samplings = prop_run_multi(dummy, 0.55, 16; PASSVALUE=passvals)
    @test size(stack) == (16, 16, 3)
    @test length(samplings) == 3
end

@testset "Phase 3 FITS/map basics" begin
    wf = prop_begin(1.0, 500e-9, 16)
    dmap = randn(16, 16)

    mktempdir() do d
        f = joinpath(d, "map.fits")
        prop_writemap(dmap, f; SAMPLING=wf.sampling_m)
        m2 = prop_readmap(wf, f; SAMPLING=wf.sampling_m)
        @test size(m2) == size(wf.field)

        f2 = joinpath(d, "img.fits")
        prop_fits_write(f2, dmap)
        r2 = prop_fits_read(f2)
        @test size(r2) == size(dmap)
    end
end
