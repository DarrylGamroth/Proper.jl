using Test
using Random

@testset "Inference and Allocation Gates" begin
    wf = prop_begin(1.0, 500e-9, 32)
    rng = MersenneTwister(123)
    dmap = rand(rng, Float32, 32, 32)

    opts = Proper.ResampleMapOptions(wf, wf.sampling_m, 16.0, 16.0)
    out = similar(dmap, Float64, size(wf.field)...)
    @test size(Proper.prop_resamplemap!(out, wf, dmap, opts)) == size(wf.field)

    # Warm-up + allocation gate for mutating resample path.
    Proper.prop_resamplemap!(out, wf, dmap, opts)
    resample_alloc = @allocated Proper.prop_resamplemap!(out, wf, dmap, opts)
    @test resample_alloc < 4_000_000

    # Scalar cubic interpolation should remain allocation-free.
    img = rand(rng, Float64, 16, 16)
    _ = prop_cubic_conv(img, 7.3, 8.1)
    @test (@allocated prop_cubic_conv(img, 7.3, 8.1)) == 0

    mktempdir() do d
        f = joinpath(d, "map.fits")
        prop_fits_write(f, dmap)

        r = @inferred AbstractArray prop_fits_read(f)
        @test size(r) == size(dmap)

        rh, hdr = @inferred Tuple{AbstractArray,Proper.FITSHeader} Proper.prop_fits_read_with_header(f)
        @test size(rh) == size(dmap)
        @test haskey(hdr, "SIMPLE")

        rm = prop_readmap(wf, f; SAMPLING=wf.sampling_m)
        @test size(rm) == size(wf.field)
    end

    dummy(λm, n; kwargs...) = prop_begin(1.0, λm, n)
    stack, samplings = prop_run_multi(dummy, 0.55, 16; PASSVALUE=[1, 2, 3])
    @test size(stack) == (16, 16, 3)
    @test eltype(samplings) <: AbstractFloat
end
