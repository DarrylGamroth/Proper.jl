using Proper
include(joinpath(@__DIR__, "multi_example.jl"))

function testmulti2(; wavelength::Real=0.6, gridsize::Integer=256, npatterns::Integer=3)
    gridsize > 0 || throw(ArgumentError("gridsize must be positive"))
    npatterns > 0 || throw(ArgumentError("npatterns must be positive"))

    nactuators = 48
    optval = [(; use_dm=true, dm=zeros(nactuators, nactuators)) for _ in 1:npatterns]
    actuator_phase = collect(range(0.0, 2pi; length=nactuators))
    x = reshape(actuator_phase, :, 1) .* ones(Float64, 1, nactuators)
    for i in 1:npatterns
        optval[i].dm .= 5e-8 .* cos.(4 .* x .* i)
    end

    fields, sampling = prop_run_multi(multi_example, wavelength, gridsize; PASSVALUE=optval)
    return fields, sampling
end

if abspath(PROGRAM_FILE) == @__FILE__
    fields, sampling = testmulti2()
    println("testmulti2: fields size = ", size(fields), ", samplings = ", sampling)
end
