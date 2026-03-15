using Base.Threads

"""Run multiple prescription passes in parallel; preserves input order."""
function prop_run_multi(routine_name, lambda0_microns, gridsize::Integer; PASSVALUE=nothing, kwargs...)
    prepared = prepare_prescription(routine_name, lambda0_microns, gridsize; PASSVALUE=PASSVALUE, kwargs...)
    return prop_run_multi(prepared)
end

"""Run a prepared prescription multiple times in parallel; preserves input order."""
function prop_run_multi(prepared::PreparedPrescription; PASSVALUE=prepared.passvalue, kwargs...)
    passes = PASSVALUE === nothing ? [nothing] : (PASSVALUE isa AbstractVector ? PASSVALUE : [PASSVALUE])
    batch = prepare_prescription_batch(prepared; pool_size=length(passes))
    return prop_run_multi(batch; PASSVALUE=passes, kwargs...)
end

"""Run a prepared prescription batch multiple times in parallel; preserves input order."""
function prop_run_multi(batch::PreparedBatch; PASSVALUE=batch.prepared.passvalue, kwargs...)
    passes = PASSVALUE === nothing ? [nothing] : (PASSVALUE isa AbstractVector ? PASSVALUE : [PASSVALUE])
    n = length(passes)
    n > 0 || throw(ArgumentError("PASSVALUE must not be empty"))
    contexts = ensure_prepared_batch_contexts!(batch, n)

    first_out, first_sampling = prop_run(batch.prepared; PASSVALUE=passes[1], context=contexts[1], kwargs...)
    first_out isa AbstractMatrix || throw(ArgumentError("prop_run_multi expects matrix outputs"))

    outT = typeof(first_out)
    sampT = typeof(float(first_sampling))
    outputs = Vector{outT}(undef, n)
    samplings = Vector{sampT}(undef, n)
    outputs[1] = first_out
    samplings[1] = float(first_sampling)
    sy, sx = size(first_out)

    @threads for i in 2:n
        out, s = prop_run(batch.prepared; PASSVALUE=passes[i], context=contexts[i], kwargs...)
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
