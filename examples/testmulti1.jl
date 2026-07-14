using Proper
using Statistics
include(joinpath(@__DIR__, "multi_example.jl"))

function testmulti1(
    ;
    lambda_min::Real=0.5,
    lambda_max::Real=0.7,
    nlambda::Integer=9,
    gridsize::Integer=256,
    npsf::Integer=256,
    final_sampling::Real=1.5e-6,
)
    nlambda > 0 || throw(ArgumentError("nlambda must be positive"))
    gridsize > 0 || throw(ArgumentError("gridsize must be positive"))
    npsf > 0 || throw(ArgumentError("npsf must be positive"))

    wavelength = collect(range(lambda_min, lambda_max; length=nlambda))

    optval = (; use_dm=true, dm=zeros(48, 48))
    optval.dm[21, 21] = 0.2e-6
    optval.dm[16, 26] = 0.2e-6

    runs = [prepare_prescription(multi_example, λ, gridsize) for λ in wavelength]
    fields, sampling = prop_run_multi(runs; PASSVALUE=fill(optval, nlambda))

    psfs = zeros(Float64, nlambda, npsf, npsf)
    for i in 1:nlambda
        mag = sampling[i] / final_sampling
        field = prop_magnify(fields[:, :, i], mag, npsf; conserve=true)
        psfs[i, :, :] .= abs2.(field)
    end

    psf = dropdims(mean(psfs; dims=1); dims=1)
    return psf
end

if abspath(PROGRAM_FILE) == @__FILE__
    psf = testmulti1()
    println("testmulti1: peak = ", maximum(psf))
end
