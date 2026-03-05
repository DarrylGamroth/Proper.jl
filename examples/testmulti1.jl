using Proper
using Statistics
include(joinpath(@__DIR__, "multi_example.jl"))

function testmulti1()
    lambda_min = 0.5
    lambda_max = 0.7
    nlambda = 9
    gridsize = 256
    npsf = 256
    final_sampling = 1.5e-6

    wavelength = collect(range(lambda_min, lambda_max; length=nlambda))

    optval = Dict("use_dm" => true, "dm" => zeros(48, 48))
    optval["dm"][21, 21] = 0.2e-6
    optval["dm"][16, 26] = 0.2e-6

    fields, sampling = prop_run_multi(multi_example, wavelength, gridsize; PASSVALUE=optval)

    psfs = zeros(Float64, nlambda, npsf, npsf)
    for i in 1:nlambda
        mag = sampling[i] / final_sampling
        field = prop_magnify(fields[:, :, i], mag, npsf)
        psfs[i, :, :] .= abs2.(field)
    end

    psf = dropdims(mean(psfs; dims=1); dims=1)
    return psf
end

if abspath(PROGRAM_FILE) == @__FILE__
    psf = testmulti1()
    println("testmulti1: peak = ", maximum(psf))
end
