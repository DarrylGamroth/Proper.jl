using Proper
include(joinpath(@__DIR__, "example_system.jl"))

function run_example(wavelength::Real, gridsize::Integer)
    state_path = joinpath(mktempdir(), "example_system.state")
    prop_init_savestate(dirname(state_path))
    psf = nothing
    sampling = 0.0
    for _ in 1:11
        psf, sampling = prop_run(example_system, wavelength * 1e6, gridsize; PASSVALUE=Dict("state_path" => state_path))
    end
    prop_end_savestate(prop_begin(1.0, wavelength, gridsize), state_path)
    return psf, sampling
end

if abspath(PROGRAM_FILE) == @__FILE__
    psf, sampling = run_example(0.5e-6, 256)
    println("run_example: sampling = ", sampling, ", peak = ", maximum(psf))
end
