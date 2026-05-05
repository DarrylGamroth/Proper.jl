@inline function _prop_run_finalize(result)
    if result isa WaveFront
        return prop_end(result)
    elseif result isa Tuple && length(result) == 2
        return result
    end
    throw(ArgumentError("Prescription must return WaveFront or (psf, sampling)"))
end

@inline _passvalue_key(key) = Symbol(lowercase(String(key)))

function _dict_to_namedtuple(dict::Dict{Symbol,Any})
    pairs_vec = collect(pairs(dict))
    keys_tuple = Tuple(first(pair) for pair in pairs_vec)
    values_tuple = Tuple(last(pair) for pair in pairs_vec)
    return NamedTuple{keys_tuple}(values_tuple)
end

function _passvalue_as_kwargs(passvalue::NamedTuple)
    dict = Dict{Symbol,Any}()
    for (key, value) in pairs(passvalue)
        kw = _passvalue_key(key)
        if haskey(dict, kw) && !isequal(dict[kw], value)
            throw(ArgumentError("PASSVALUE contains conflicting values for keyword $(kw)"))
        end
        dict[kw] = value
    end
    return _dict_to_namedtuple(dict)
end

function _passvalue_as_kwargs(passvalue::AbstractDict)
    dict = Dict{Symbol,Any}()
    for (key, value) in pairs(passvalue)
        kw = _passvalue_key(key)
        if haskey(dict, kw) && !isequal(dict[kw], value)
            throw(ArgumentError("PASSVALUE contains conflicting values for keyword $(kw)"))
        end
        dict[kw] = value
    end
    return _dict_to_namedtuple(dict)
end

function _merge_passvalue_kwargs(passvalue_kwargs::NamedTuple, kwargs::NamedTuple)
    dict = Dict{Symbol,Any}()
    for (key, value) in pairs(passvalue_kwargs)
        dict[key] = value
    end
    for (key, value) in pairs(kwargs)
        if haskey(dict, key) && !isequal(dict[key], value)
            throw(ArgumentError("PASSVALUE and explicit keyword arguments both define $(key) with different values"))
        end
        dict[key] = value
    end
    return _dict_to_namedtuple(dict)
end

@inline function _call_prescription(fn::F, λm, gridsize::Integer, ::Nothing, kwargs::NamedTuple) where {F<:Function}
    return fn(λm, gridsize; kwargs...)
end

@inline function _call_prescription(fn::F, λm, gridsize::Integer, passvalue::Union{NamedTuple,AbstractDict}, kwargs::NamedTuple) where {F<:Function}
    native_kwargs = _merge_passvalue_kwargs(_passvalue_as_kwargs(passvalue), kwargs)
    return fn(λm, gridsize; native_kwargs...)
end

@inline function _call_prescription(fn::F, λm, gridsize::Integer, passvalue, kwargs::NamedTuple) where {F<:Function}
    try
        return fn(λm, gridsize, passvalue; kwargs...)
    catch err
        if err isa MethodError && err.f === fn
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

@inline function _call_hot_prescription(fn::F, λm, gridsize::Integer, kwargs::NamedTuple) where {F<:Function}
    return fn(λm, gridsize; kwargs...)
end

@inline function _prop_run_hot(call::PreparedHotCall{F,T,CTX,KW,HotCallContextActive}) where {F,T,CTX,KW}
    result = with_run_context(call.context) do
        _call_hot_prescription(call.routine, call.wavelength_m, call.gridsize, call.kwargs)
    end
    return _prop_run_finalize(result)
end

@inline function _prop_run_hot(call::PreparedHotCall{F,T,CTX,KW,HotCallContextInactive}) where {F,T,CTX,KW}
    result = _call_hot_prescription(call.routine, call.wavelength_m, call.gridsize, call.kwargs)
    return _prop_run_finalize(result)
end

"""
    prop_run_hot(call::PreparedHotCall)

Execute a pre-bound native Julia prescription hot call.

This bypasses per-call model asset resolution, slot lookup, and keyword merging.
Construct `call` with `prepare_hot_call`.
"""
function prop_run_hot(call::PreparedHotCall)
    return _prop_run_hot(call)
end

"""
    prop_run(routine_name, lambda0_microns, gridsize; PASSVALUE=nothing, context=nothing, kwargs...)
    prop_run(prepared::PreparedPrescription; PASSVALUE=prepared.passvalue, context=prepared.context, kwargs...)
    prop_run(batch::PreparedBatch; PASSVALUE=batch.prepared.passvalue, slot=1, kwargs...)
    prop_run(model::PreparedModel; PASSVALUE=model.prepared.passvalue, slot=1, kwargs...)

Execute a PROPER prescription and return `(psf, pixscale)`.

# Arguments
- `routine_name`: function object or global name of the prescription
- `lambda0_microns`: wavelength in microns, or the normalized wavelength
  already stored in a prepared object
- `gridsize`: computational grid dimension; should be a power of two for the
  usual FFT-based paths

# Keywords
- `PASSVALUE`: compatibility adapter for upstream-style optional inputs.
  `NamedTuple` and dictionary values are normalized into native Julia keyword
  arguments before the prescription is called; other values are forwarded using
  the legacy positional/`PASSVALUE` calling style.
- `context`: explicit `RunContext` to reuse workspace and backend state
- additional keyword arguments are passed through to the prescription

# Returns
- `(psf, pixscale)`, where `psf` is the prescription output and `pixscale` is
  the output sampling in meters per pixel.

# Notes
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

prop_run(call::PreparedHotCall) = prop_run_hot(call)
