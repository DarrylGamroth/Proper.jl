using Base.Threads

@inline function _run_single_pass(routine_name, lambda0_microns, gridsize::Integer, pass, kwargs)
    return pass === nothing ?
        prop_run(routine_name, lambda0_microns, gridsize; kwargs...) :
        prop_run(routine_name, lambda0_microns, gridsize; PASSVALUE=pass, kwargs...)
end

"""Run multiple prescription passes in parallel; preserves input order."""
function prop_run_multi(routine_name, lambda0_microns, gridsize::Integer; PASSVALUE=nothing, kwargs...)
    passes = PASSVALUE === nothing ? [nothing] : (PASSVALUE isa AbstractVector ? PASSVALUE : [PASSVALUE])
    n = length(passes)
    n > 0 || throw(ArgumentError("PASSVALUE must not be empty"))

    first_out, first_sampling = _run_single_pass(routine_name, lambda0_microns, gridsize, passes[1], kwargs)
    first_out isa AbstractMatrix || throw(ArgumentError("prop_run_multi expects matrix outputs"))

    outT = typeof(first_out)
    sampT = typeof(float(first_sampling))
    outputs = Vector{outT}(undef, n)
    samplings = Vector{sampT}(undef, n)
    outputs[1] = first_out
    samplings[1] = float(first_sampling)
    sy, sx = size(first_out)

    @threads for i in 2:n
        out, s = _run_single_pass(routine_name, lambda0_microns, gridsize, passes[i], kwargs)
        out isa outT || throw(ArgumentError("All prop_run_multi outputs must have the same type"))
        size(out) == (sy, sx) || throw(ArgumentError("All prop_run_multi outputs must have the same size"))
        outputs[i] = out
        samplings[i] = float(s)
    end

    stack = Array{eltype(first_out),3}(undef, sy, sx, n)
    @inbounds for i in 1:n
        stack[:, :, i] = outputs[i]
    end
    return stack, samplings
end
