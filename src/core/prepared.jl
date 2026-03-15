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
    kwargs...,
)
    prepared = prepare_prescription(routine_name, lambda0_microns, gridsize; kwargs...)
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

function reset_prepared_batch!(batch::PreparedBatch)
    for ctx in batch.contexts
        ctx === nothing && continue
        reset_workspace!(ctx.workspace)
    end
    return batch
end

function prepare_model(
    routine_name,
    lambda0_microns::Real,
    gridsize::Integer;
    name=routine_name,
    assets=nothing,
    pool_size::Integer=Base.Threads.nthreads(),
    kwargs...,
)
    prepared = prepare_prescription(routine_name, lambda0_microns, gridsize; kwargs...)
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
    return model
end

function prepare_prescription(
    routine_name,
    lambda0_microns::Real,
    gridsize::Integer;
    PASSVALUE=nothing,
    context::Union{Nothing,RunContext}=nothing,
    kwargs...,
)
    λm = float(lambda0_microns) * 1e-6
    fn = resolve_prescription_routine(routine_name)
    return PreparedPrescription(fn, λm, Int(gridsize), context, (; kwargs...), PASSVALUE)
end
