using Proper
include(joinpath(@__DIR__, "multi_example.jl"))

function testmulti2()
    wavelength = 0.6
    gridsize = 256

    npatterns = 3
    optval = [Dict("use_dm" => true, "dm" => zeros(gridsize, gridsize)) for _ in 1:npatterns]
    x = (collect(0:(gridsize - 1)) ./ (gridsize - 1) .* (2pi)) .* ones(Float64, 1, gridsize)
    for i in 1:npatterns
        optval[i]["dm"] .= 5e-8 .* cos.(4 .* x .* i)
    end

    fields, sampling = prop_run_multi(multi_example, wavelength, gridsize; PASSVALUE=optval)
    return fields, sampling
end

if abspath(PROGRAM_FILE) == @__FILE__
    fields, sampling = testmulti2()
    println("testmulti2: fields size = ", size(fields), ", samplings = ", sampling)
end
