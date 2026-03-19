@inline resolve_prescription_routine(routine::F) where {F<:Function} = routine
@inline resolve_prescription_routine(routine::Symbol) = getfield(Main, routine)
@inline resolve_prescription_routine(routine::AbstractString) = getfield(Main, Symbol(routine))

function resolve_prescription_routine(routine)
    throw(ArgumentError("routine_name must be Function, Symbol, or String"))
end

struct PreparedPrescription{F,T<:AbstractFloat,CTX,KW,P}
    routine::F
    wavelength_m::T
    gridsize::Int
    context::CTX
    kwargs::KW
    passvalue::P
end

mutable struct PreparedBatch{P,CV<:AbstractVector}
    prepared::P
    contexts::CV
end

mutable struct PreparedAssetPool{F,CV<:AbstractVector}
    factory::F
    cache::CV
end

mutable struct PreparedModel{ID,P,B,A}
    name::ID
    prepared::P
    batch::B
    assets::A
end

@inline prepared_context(prepared::PreparedPrescription) = prepared.context
@inline prepared_prescription(batch::PreparedBatch) = batch.prepared
@inline prepared_prescription(model::PreparedModel) = model.prepared
@inline prepared_batch(model::PreparedModel) = model.batch
@inline prepared_assets(model::PreparedModel) = model.assets

@inline function _prepared_context_with_precision(
    context::Union{Nothing,RunContext},
    ::Nothing,
)
    return context
end

@inline function _prepared_context_with_precision(
    ::Nothing,
    precision::Type{T},
) where {T<:AbstractFloat}
    return RunContext(Matrix{T})
end

@inline function _prepared_context_with_precision(
    context::RunContext,
    precision::Type{T},
) where {T<:AbstractFloat}
    ctxT = workspace_float_type(context.workspace)
    ctxT === T || throw(ArgumentError("context workspace precision $(ctxT) does not match requested precision $(T)"))
    return context
end

function prepare_asset_pool(factory; pool_size::Integer=Base.Threads.nthreads())
    pool_size > 0 || throw(ArgumentError("pool_size must be positive"))
    cache = Vector{Any}(undef, pool_size)
    fill!(cache, nothing)
    return PreparedAssetPool(factory, cache)
end

function prepared_contexts(prepared::PreparedPrescription, n::Integer)
    n >= 0 || throw(ArgumentError("n must be non-negative"))
    ctx = prepared.context
    ctx === nothing && return fill(nothing, n)

    contexts = Vector{RunContext}(undef, n)
    for i in 1:n
        contexts[i] = fresh_context(ctx; rng=fork_rng(ctx.rng, i))
    end
    return contexts
end

function prepare_prescription_batch(prepared::PreparedPrescription; pool_size::Integer=Base.Threads.nthreads())
    pool_size > 0 || throw(ArgumentError("pool_size must be positive"))
    return PreparedBatch(prepared, prepared_contexts(prepared, pool_size))
end

function prepare_prescription_batch(
    routine_name,
    lambda0_microns::Real,
    gridsize::Integer;
    pool_size::Integer=Base.Threads.nthreads(),
    precision::Union{Nothing,Type{<:AbstractFloat}}=nothing,
    kwargs...,
)
    prepared = prepare_prescription(routine_name, lambda0_microns, gridsize; precision=precision, kwargs...)
    return prepare_prescription_batch(prepared; pool_size=pool_size)
end

function ensure_prepared_batch_contexts!(batch::PreparedBatch, n::Integer)
    n >= 0 || throw(ArgumentError("n must be non-negative"))
    current = length(batch.contexts)
    current >= n && return batch.contexts

    resize!(batch.contexts, n)
    ctx = prepared_context(batch.prepared)
    if ctx === nothing
        @inbounds for i in current + 1:n
            batch.contexts[i] = nothing
        end
    else
        @inbounds for i in current + 1:n
            batch.contexts[i] = fresh_context(ctx; rng=fork_rng(ctx.rng, i))
        end
    end
    return batch.contexts
end

function ensure_prepared_asset_pool!(pool::PreparedAssetPool, n::Integer)
    n >= 0 || throw(ArgumentError("n must be non-negative"))
    current = length(pool.cache)
    current >= n && return pool.cache

    resize!(pool.cache, n)
    @inbounds for i in current + 1:n
        pool.cache[i] = nothing
    end
    return pool.cache
end

function reset_prepared_batch!(batch::PreparedBatch)
    for ctx in batch.contexts
        ctx === nothing && continue
        reset_workspace!(ctx.workspace)
    end
    return batch
end

function reset_prepared_assets!(pool::PreparedAssetPool)
    fill!(pool.cache, nothing)
    return pool
end

@inline reset_prepared_assets!(::Nothing) = nothing
@inline reset_prepared_assets!(assets) = assets

@inline function _invoke_asset_factory(factory, slot::Integer, model::PreparedModel)
    if applicable(factory, slot, model)
        return factory(slot, model)
    elseif applicable(factory, slot)
        return factory(slot)
    elseif applicable(factory, model)
        return factory(model)
    end
    return factory()
end

@inline _resolve_model_assets(::Nothing, ::PreparedModel, ::Integer) = nothing
@inline _resolve_model_assets(assets, ::PreparedModel, ::Integer) = assets

function _resolve_model_assets(pool::PreparedAssetPool, model::PreparedModel, slot::Integer)
    slot > 0 || throw(ArgumentError("slot must be positive"))
    cache = ensure_prepared_asset_pool!(pool, slot)
    asset = cache[slot]
    if asset === nothing
        asset = _invoke_asset_factory(pool.factory, slot, model)
        cache[slot] = asset
    end
    return asset
end

@inline prepared_assets(model::PreparedModel, slot::Integer) = _resolve_model_assets(model.assets, model, slot)

function prepare_model(
    routine_name,
    lambda0_microns::Real,
    gridsize::Integer;
    name=routine_name,
    assets=nothing,
    pool_size::Integer=Base.Threads.nthreads(),
    precision::Union{Nothing,Type{<:AbstractFloat}}=nothing,
    kwargs...,
)
    prepared = prepare_prescription(routine_name, lambda0_microns, gridsize; precision=precision, kwargs...)
    batch = prepare_prescription_batch(prepared; pool_size=pool_size)
    return PreparedModel(name, prepared, batch, assets)
end

function prepare_model(
    name,
    routine_name,
    lambda0_microns::Real,
    gridsize::Integer;
    assets=nothing,
    pool_size::Integer=Base.Threads.nthreads(),
    precision::Union{Nothing,Type{<:AbstractFloat}}=nothing,
    kwargs...,
)
    prepared = prepare_prescription(routine_name, lambda0_microns, gridsize; precision=precision, kwargs...)
    batch = prepare_prescription_batch(prepared; pool_size=pool_size)
    return PreparedModel(name, prepared, batch, assets)
end

function prepare_model(
    prepared::PreparedPrescription;
    name=prepared.routine,
    assets=nothing,
    pool_size::Integer=Base.Threads.nthreads(),
)
    batch = prepare_prescription_batch(prepared; pool_size=pool_size)
    return PreparedModel(name, prepared, batch, assets)
end

function reset_prepared_model!(model::PreparedModel)
    reset_prepared_batch!(model.batch)
    reset_prepared_assets!(model.assets)
    return model
end

function prepare_prescription(
    routine_name,
    lambda0_microns::Real,
    gridsize::Integer;
    PASSVALUE=nothing,
    context::Union{Nothing,RunContext}=nothing,
    precision::Union{Nothing,Type{<:AbstractFloat}}=nothing,
    kwargs...,
)
    WT = isnothing(precision) ? float(typeof(lambda0_microns)) : precision
    λm = WT(lambda0_microns) * WT(1e-6)
    fn = resolve_prescription_routine(routine_name)
    ctx = _prepared_context_with_precision(context, precision)
    return PreparedPrescription(fn, λm, Int(gridsize), ctx, (; kwargs...), PASSVALUE)
end
