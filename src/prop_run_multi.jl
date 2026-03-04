using Base.Threads

"""Run multiple prescription passes in parallel; preserves input order."""
function prop_run_multi(routine_name, lambda0_microns, gridsize::Integer; PASSVALUE=nothing, kwargs...)
    passes = PASSVALUE === nothing ? [nothing] : (PASSVALUE isa AbstractVector ? PASSVALUE : [PASSVALUE])
    n = length(passes)
    outputs = Vector{Any}(undef, n)
    samplings = Vector{Float64}(undef, n)

    @threads for i in 1:n
        out, s = passes[i] === nothing ?
            prop_run(routine_name, lambda0_microns, gridsize; kwargs...) :
            prop_run(routine_name, lambda0_microns, gridsize; PASSVALUE=passes[i], kwargs...)
        outputs[i] = out
        samplings[i] = float(s)
    end

    first_out = outputs[1]
    stack = Array{eltype(first_out),3}(undef, size(first_out, 1), size(first_out, 2), n)
    @inbounds for i in 1:n
        stack[:, :, i] = outputs[i]
    end
    return stack, samplings
end
