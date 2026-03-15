@inline function _prop_run_finalize(result)
    if result isa WaveFront
        return prop_end(result)
    elseif result isa Tuple && length(result) == 2
        return result
    end
    throw(ArgumentError("Prescription must return WaveFront or (psf, sampling)"))
end

@inline function _call_prescription(fn::F, λm, gridsize::Integer, ::Nothing, kwargs::NamedTuple) where {F<:Function}
    return fn(λm, gridsize; kwargs...)
end

@inline function _call_prescription(fn::F, λm, gridsize::Integer, passvalue, kwargs::NamedTuple) where {F<:Function}
    try
        return fn(λm, gridsize, passvalue; kwargs...)
    catch err
        if err isa MethodError
            return fn(λm, gridsize; PASSVALUE=passvalue, kwargs...)
        end
        rethrow(err)
    end
end

function _prop_run_resolved(
    fn::F,
    λm,
    gridsize::Integer,
    passvalue,
    context::Union{Nothing,RunContext},
    kwargs::NamedTuple,
) where {F<:Function}
    result = with_run_context(context) do
        _call_prescription(fn, λm, gridsize, passvalue, kwargs)
    end
    return _prop_run_finalize(result)
end

"""Run a prescription function by function object or global name."""
function prop_run(
    routine_name,
    lambda0_microns::Real,
    gridsize::Integer;
    PASSVALUE=nothing,
    context::Union{Nothing,RunContext}=nothing,
    kwargs...,
)
    prepared = prepare_prescription(
        routine_name,
        lambda0_microns,
        gridsize;
        PASSVALUE=PASSVALUE,
        context=context,
        kwargs...,
    )
    return prop_run(prepared)
end

function prop_run(
    prepared::PreparedPrescription;
    PASSVALUE=prepared.passvalue,
    context::Union{Nothing,RunContext}=prepared.context,
    kwargs...,
)
    merged_kwargs = merge(prepared.kwargs, (; kwargs...))
    return _prop_run_resolved(prepared.routine, prepared.wavelength_m, prepared.gridsize, PASSVALUE, context, merged_kwargs)
end

function prop_run(
    batch::PreparedBatch;
    PASSVALUE=batch.prepared.passvalue,
    slot::Integer=1,
    kwargs...,
)
    slot > 0 || throw(ArgumentError("slot must be positive"))
    contexts = ensure_prepared_batch_contexts!(batch, slot)
    return prop_run(batch.prepared; PASSVALUE=PASSVALUE, context=contexts[slot], kwargs...)
end
