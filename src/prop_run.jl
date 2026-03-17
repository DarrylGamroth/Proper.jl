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

"""
    prop_run(routine_name, lambda0_microns, gridsize; PASSVALUE=nothing, context=nothing, kwargs...)
    prop_run(prepared::PreparedPrescription; PASSVALUE=prepared.passvalue, context=prepared.context, kwargs...)
    prop_run(batch::PreparedBatch; PASSVALUE=batch.prepared.passvalue, slot=1, kwargs...)
    prop_run(model::PreparedModel; PASSVALUE=model.prepared.passvalue, slot=1, kwargs...)

Execute a PROPER prescription and return `(psf, pixscale)`.

Outputs:
- `psf`: image or field returned by the prescription, after `prop_end` if the
  prescription returns a `WaveFront`
- `pixscale`: output sampling in meters per pixel

Required inputs:
- `routine_name`: function object or global name of the prescription
- `lambda0_microns`: wavelength in microns, or the normalized wavelength
  already stored in a prepared object
- `gridsize`: computational grid dimension; should be a power of two for the
  usual FFT-based paths

Optional inputs:
- `PASSVALUE`: value forwarded to the prescription in the familiar PROPER
  calling style
- `context`: explicit `RunContext` to reuse workspace and backend state
- additional keyword arguments are passed through to the prescription

Notes:
- Prepared forms reuse normalized arguments, run contexts, and optional
  prepared assets while preserving the public return contract.
- Prescriptions may return either a `WaveFront` or `(psf, sampling)`.
"""
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

function prop_run(
    model::PreparedModel;
    PASSVALUE=model.prepared.passvalue,
    slot::Integer=1,
    kwargs...,
)
    assets = prepared_assets(model, slot)
    merged_kwargs = assets === nothing ? (; kwargs...) :
        (assets isa NamedTuple ? merge(assets, (; kwargs...)) : merge((; assets=assets), (; kwargs...)))
    return prop_run(model.batch; PASSVALUE=PASSVALUE, slot=slot, merged_kwargs...)
end
