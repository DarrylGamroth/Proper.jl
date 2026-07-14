using Base.Threads

abstract type MultiRunExecStyle end
struct MultiRunThreadedExecStyle <: MultiRunExecStyle end
struct MultiRunSerialExecStyle <: MultiRunExecStyle end

const MultiRunItem = Union{PreparedPrescription,PreparedBatch,PreparedModel,PreparedRun}

@inline _multi_run_exec_style(::CPUBackend) = MultiRunThreadedExecStyle()
@inline _multi_run_exec_style(::BackendStyle) = MultiRunSerialExecStyle()
@inline _multi_run_exec_style(::Nothing) = MultiRunThreadedExecStyle()
@inline _multi_run_exec_style(ctx::RunContext) = _multi_run_exec_style(ctx.backend)

@inline _multi_run_context(prepared::PreparedPrescription) = prepared.context
@inline _multi_run_context(batch::PreparedBatch) = prepared_context(batch, 1)
@inline _multi_run_context(model::PreparedModel) = prepared_context(model, 1)
@inline _multi_run_context(prepared::PreparedRun) = prepared.context

function _multi_run_exec_style(contexts::AbstractVector)
    @inbounds for i in eachindex(contexts)
        ctx = contexts[i]
        _multi_run_exec_style(ctx) isa MultiRunSerialExecStyle && return MultiRunSerialExecStyle()
        if ctx isa RunContext
            for j in firstindex(contexts):(i - 1)
                previous = contexts[j]
                previous isa RunContext && previous.workspace === ctx.workspace &&
                    return MultiRunSerialExecStyle()
            end
        end
    end
    return MultiRunThreadedExecStyle()
end

function _multi_run_exec_style(
    runs::AbstractVector{<:MultiRunItem},
)
    @inbounds for i in eachindex(runs)
        ctx = _multi_run_context(runs[i])
        _multi_run_exec_style(ctx) isa MultiRunSerialExecStyle && return MultiRunSerialExecStyle()
        if ctx isa RunContext
            for j in firstindex(runs):(i - 1)
                previous = _multi_run_context(runs[j])
                previous isa RunContext && previous.workspace === ctx.workspace &&
                    return MultiRunSerialExecStyle()
            end
        end
    end
    return MultiRunThreadedExecStyle()
end

struct MultiRunRepeatedPass{P}
    value::P
    length::Int
end

Base.length(passes::MultiRunRepeatedPass) = passes.length

@inline function Base.getindex(passes::MultiRunRepeatedPass, i::Integer)
    @boundscheck checkbounds(Base.OneTo(passes.length), i)
    return passes.value
end

@inline _multi_run_vector_passes(passvalue::AbstractVector, ::Integer) = passvalue
@inline _multi_run_vector_passes(passvalue, n::Integer) = MultiRunRepeatedPass(passvalue, Int(n))

struct MultiRunItemRunner{KW}
    kwargs::KW
end

@inline function (runner::MultiRunItemRunner)(run::PreparedRun, pass, ::Integer)
    pass === nothing || throw(ArgumentError("Prepared runs do not accept PASSVALUE; bind native Julia keywords with prepare_run"))
    isempty(runner.kwargs) || throw(ArgumentError("Prepared runs do not accept execution keywords; bind them with prepare_run"))
    return prop_run(run)
end
@inline function (runner::MultiRunItemRunner)(
    run::Union{PreparedPrescription,PreparedBatch,PreparedModel},
    pass,
    ::Integer,
)
    return prop_run(run; PASSVALUE=pass, runner.kwargs...)
end

@inline _multi_passes(passvalue) = passvalue === nothing ? [nothing] : (passvalue isa AbstractVector ? passvalue : [passvalue])

@inline function _allocate_output_stack(first_out::AbstractMatrix, sy::Integer, sx::Integer, n::Integer)
    return similar(first_out, eltype(first_out), sy, sx, n)
end

@inline function _store_output_slice!(stack, i::Integer, out::AbstractMatrix)
    copyto!(selectdim(stack, 3, i), out)
    return stack
end

function _finish_multi_run_serial(
    first_out::O,
    first_sampling,
    items::AbstractVector,
    passes,
    runner::F,
) where {O<:AbstractMatrix,F}
    n = length(items)
    sampT = typeof(float(first_sampling))
    sy, sx = size(first_out)
    stack = _allocate_output_stack(first_out, sy, sx, n)
    samplings = Vector{sampT}(undef, n)
    _store_output_slice!(stack, 1, first_out)
    samplings[1] = float(first_sampling)

    @inbounds for i in 2:n
        out, s = runner(items[i], passes[i], i)
        out isa O || throw(ArgumentError("All prop_run_multi outputs must have the same type"))
        size(out) == (sy, sx) || throw(ArgumentError("All prop_run_multi outputs must have the same size"))
        _store_output_slice!(stack, i, out)
        samplings[i] = float(s)
    end
    return stack, samplings
end

function _prop_run_multi_items(
    ::MultiRunSerialExecStyle,
    items::AbstractVector,
    passes,
    runner::F,
) where {F}
    n = length(items)
    n == length(passes) || throw(ArgumentError("PASSVALUE length must match number of prepared runs"))
    n > 0 || throw(ArgumentError("prepared run collection must not be empty"))

    first_out, first_sampling = runner(items[1], passes[1], 1)
    first_out isa AbstractMatrix || throw(ArgumentError("prop_run_multi expects matrix outputs"))
    return _finish_multi_run_serial(first_out, first_sampling, items, passes, runner)
end

@inline function _validate_multi_run_destinations(
    stack::AbstractArray{<:Number,3},
    samplings::AbstractVector{<:AbstractFloat},
    n::Integer,
)
    n > 0 || throw(ArgumentError("prepared run collection must not be empty"))
    size(stack, 3) == n || throw(ArgumentError(
        "output stack depth $(size(stack, 3)) must match number of runs $(n)",
    ))
    length(samplings) == n || throw(ArgumentError(
        "sampling vector length $(length(samplings)) must match number of runs $(n)",
    ))
    backend_style(typeof(samplings)) isa CPUBackend || throw(ArgumentError(
        "sampling metadata must use a host vector",
    ))
    return Int(n)
end

@inline _is_exact_output_slice(::AbstractMatrix, ::AbstractArray, ::Integer) = false

@inline function _is_exact_output_slice(out::SubArray, stack::AbstractArray, i::Integer)
    parent(out) === stack || return false
    indices = parentindices(out)
    length(indices) == 3 || return false
    return indices[1] == Base.Slice(axes(stack, 1)) &&
        indices[2] == Base.Slice(axes(stack, 2)) &&
        indices[3] == i
end

@inline function _store_output_slice!(
    stack::AbstractArray{<:Number,3},
    samplings::AbstractVector{<:AbstractFloat},
    i::Integer,
    out,
    sampling,
)
    out isa AbstractMatrix || throw(ArgumentError("prop_run_multi! expects matrix outputs"))
    size(out) == (size(stack, 1), size(stack, 2)) || throw(ArgumentError(
        "prop_run_multi! output size $(size(out)) does not match destination size $(size(stack)[1:2])",
    ))
    eltype(out) === eltype(stack) || throw(ArgumentError(
        "prop_run_multi! output eltype $(eltype(out)) does not match destination eltype $(eltype(stack))",
    ))
    same_backend_style(typeof(stack), typeof(out)) || throw(ArgumentError(
        "prop_run_multi! output and destination must use the same backend",
    ))

    _is_exact_output_slice(out, stack, i) || copyto!(selectdim(stack, 3, i), out)
    samplings[i] = float(sampling)
    return stack
end

@inline _multi_run_mutating_exec_style(
    style::MultiRunSerialExecStyle,
    ::AbstractArray,
    ::AbstractVector,
) = style

@inline _multi_run_mutating_exec_style(
    style::MultiRunThreadedExecStyle,
    ::StridedArray,
    ::StridedVector,
) = style

@inline _multi_run_mutating_exec_style(
    ::MultiRunThreadedExecStyle,
    ::AbstractArray,
    ::AbstractVector,
) = MultiRunSerialExecStyle()

function _prop_run_multi_items!(
    ::MultiRunSerialExecStyle,
    stack::AbstractArray{<:Number,3},
    samplings::AbstractVector{<:AbstractFloat},
    items::AbstractVector,
    passes,
    runner::F,
) where {F}
    n = _validate_multi_run_destinations(stack, samplings, length(items))
    n == length(passes) || throw(ArgumentError("PASSVALUE length must match number of prepared runs"))

    @inbounds for i in 1:n
        out, sampling = runner(items[i], passes[i], i)
        _store_output_slice!(stack, samplings, i, out, sampling)
    end
    return stack, samplings
end

function _prop_run_multi_items!(
    ::MultiRunThreadedExecStyle,
    stack::AbstractArray{<:Number,3},
    samplings::AbstractVector{<:AbstractFloat},
    items::AbstractVector,
    passes,
    runner::F,
) where {F}
    n = _validate_multi_run_destinations(stack, samplings, length(items))
    n == length(passes) || throw(ArgumentError("PASSVALUE length must match number of prepared runs"))
    if Threads.nthreads() == 1 || n == 1
        return _prop_run_multi_items!(
            MultiRunSerialExecStyle(),
            stack,
            samplings,
            items,
            passes,
            runner,
        )
    end

    inner_min_elems = n >= Threads.nthreads() ? 262_144 : 0
    inner_policy = CPUInnerKernelPolicy(true, inner_min_elems)
    @threads for i in 1:n
        out, sampling = with_cpu_inner_kernel_policy(inner_policy) do
            runner(items[i], passes[i], i)
        end
        _store_output_slice!(stack, samplings, i, out, sampling)
    end
    return stack, samplings
end

struct MultiRunPreparedRunner{P,KW}
    prepared::P
    kwargs::KW
end

@inline function (runner::MultiRunPreparedRunner)(ctx, pass, ::Integer)
    return prop_run(
        runner.prepared;
        PASSVALUE=pass,
        context=ctx,
        runner.kwargs...,
    )
end

struct MultiRunModelRunner{M,KW}
    model::M
    kwargs::KW
end

@inline function (runner::MultiRunModelRunner)(::Integer, pass, slot::Integer)
    return prop_run(runner.model; PASSVALUE=pass, slot=slot, runner.kwargs...)
end

function _finish_multi_run_threaded(
    first_out::O,
    first_sampling,
    items::AbstractVector,
    passes,
    runner::F,
    inner_policy::CPUInnerKernelPolicy,
) where {O<:AbstractMatrix,F}
    n = length(items)
    sampT = typeof(float(first_sampling))
    sy, sx = size(first_out)
    results = Vector{O}(undef, n)
    samplings = Vector{sampT}(undef, n)
    results[1] = first_out
    samplings[1] = float(first_sampling)

    @threads for i in 2:n
        out, s = with_cpu_inner_kernel_policy(inner_policy) do
            runner(items[i], passes[i], i)
        end
        out isa O || throw(ArgumentError("All prop_run_multi outputs must have the same type"))
        size(out) == (sy, sx) || throw(ArgumentError("All prop_run_multi outputs must have the same size"))
        results[i] = out
        samplings[i] = float(s)
    end

    stack = _allocate_output_stack(first_out, sy, sx, n)
    @inbounds for i in eachindex(results)
        _store_output_slice!(stack, i, results[i])
    end

    return stack, samplings
end

function _prop_run_multi_items(
    ::MultiRunThreadedExecStyle,
    items::AbstractVector,
    passes,
    runner::F,
) where {F}
    n = length(items)
    n == length(passes) || throw(ArgumentError("PASSVALUE length must match number of prepared runs"))
    n > 0 || throw(ArgumentError("prepared run collection must not be empty"))

    if Threads.nthreads() == 1 || n == 1
        return _prop_run_multi_items(MultiRunSerialExecStyle(), items, passes, runner)
    end

    inner_min_elems = n >= Threads.nthreads() ? 262_144 : 0
    inner_policy = CPUInnerKernelPolicy(true, inner_min_elems)
    first_out, first_sampling = with_cpu_inner_kernel_policy(inner_policy) do
        runner(items[1], passes[1], 1)
    end
    first_out isa AbstractMatrix || throw(ArgumentError("prop_run_multi expects matrix outputs"))
    return _finish_multi_run_threaded(first_out, first_sampling, items, passes, runner, inner_policy)
end

"""
    prop_run_multi(routine_name, lambda0_microns, gridsize; PASSVALUE=nothing, kwargs...)
    prop_run_multi(prepared::PreparedPrescription; PASSVALUE=prepared.passvalue, kwargs...)
    prop_run_multi(batch::PreparedBatch; PASSVALUE=batch.prepared.passvalue, kwargs...)
    prop_run_multi(model::PreparedModel; PASSVALUE=model.prepared.passvalue, kwargs...)
    prop_run_multi(runs::AbstractVector{<:MultiRunItem}; PASSVALUE=nothing, kwargs...)

Execute multiple prescription instances and preserve input order.

# Arguments
- `routine_name`: function object or global name of the prescription
- `lambda0_microns`: scalar or prepared wavelength; repeated runs are driven by
  `PASSVALUE`
- `gridsize`: computational grid size shared by all runs

# Keywords
- `PASSVALUE`: scalar or vector of values to forward to the prescription; when
  a vector is supplied, each entry is run independently
- additional keyword arguments are passed through to each prescription call

# Returns
- `(stack, samplings)`, where `stack[:, :, i]` is the `i`th output and
  `samplings[i]` is its sampling in meters per pixel.

# Notes
- Prepared batch and model forms reuse per-slot run contexts so repeated
  multi-run workloads avoid rebuilding core workspace state.
- Vector forms support wavelength sweeps or mixed prepared runs while
  preserving backend-native stacked outputs where feasible.
- CPU runs use Julia threads when multiple independent contexts are available.
  Accelerator and unknown backends use ordered serial host submission by
  default, avoiding contention between host tasks that target one device.
- Output order is deterministic and matches the input order.
"""
function prop_run_multi(routine_name, lambda0_microns, gridsize::Integer; PASSVALUE=nothing, kwargs...)
    prepared = prepare_prescription(routine_name, lambda0_microns, gridsize; PASSVALUE=PASSVALUE, kwargs...)
    return prop_run_multi(prepared)
end

"""Run a prepared prescription multiple times; preserves input order."""
function prop_run_multi(prepared::PreparedPrescription; PASSVALUE=prepared.passvalue, kwargs...)
    passes = _multi_passes(PASSVALUE)
    batch = prepare_prescription_batch(prepared; pool_size=length(passes))
    return prop_run_multi(batch; PASSVALUE=passes, kwargs...)
end

"""Run a prepared prescription batch multiple times; preserves input order."""
function prop_run_multi(batch::PreparedBatch; PASSVALUE=batch.prepared.passvalue, kwargs...)
    passes = _multi_passes(PASSVALUE)
    n = length(passes)
    n > 0 || throw(ArgumentError("PASSVALUE must not be empty"))
    contexts = ensure_prepared_batch_contexts!(batch, n)
    items = @view contexts[1:n]
    style = _multi_run_exec_style(items)
    return _prop_run_multi_items(style, items, passes, (ctx, pass, _) -> begin
        prop_run(batch.prepared; PASSVALUE=pass, context=ctx, kwargs...)
    end)
end

function prop_run_multi(model::PreparedModel; PASSVALUE=model.prepared.passvalue, kwargs...)
    passes = _multi_passes(PASSVALUE)
    n = length(passes)
    n > 0 || throw(ArgumentError("PASSVALUE must not be empty"))
    contexts = ensure_prepared_batch_contexts!(model.batch, n)
    model.assets isa PreparedAssetPool && ensure_prepared_asset_pool!(model.assets, n)
    style = _multi_run_exec_style(@view contexts[1:n])
    return _prop_run_multi_items(style, Base.OneTo(n), passes, (_, pass, slot) -> begin
        prop_run(model; PASSVALUE=pass, slot=slot, kwargs...)
    end)
end

function prop_run_multi(
    runs::AbstractVector{<:MultiRunItem};
    PASSVALUE=nothing,
    kwargs...,
)
    passes = _multi_run_vector_passes(PASSVALUE, length(runs))
    style = _multi_run_exec_style(runs)
    return _prop_run_multi_items(style, runs, passes, MultiRunItemRunner((; kwargs...)))
end

"""
    prop_run_multi!(stack, samplings, routine_name, lambda0_microns, gridsize; PASSVALUE=nothing, kwargs...)
    prop_run_multi!(stack, samplings, prepared; PASSVALUE=prepared.passvalue, kwargs...)
    prop_run_multi!(stack, samplings, runs; PASSVALUE=nothing, kwargs...)

Execute a multi-run workload into caller-owned output and sampling storage.

`stack` must be a three-dimensional numeric array whose third-axis length
matches the run count. `samplings` must be a host floating-point vector of the
same length. Every prescription output must match the stack's first two
dimensions, element type, and backend exactly.

Dense strided CPU destinations may be populated concurrently when run contexts
own independent workspaces. Packed or custom destinations use serial execution
because logically separate slices may share physical storage words.

For the lowest steady-state allocation overhead, bind independent wavefront and
output buffers with `prepare_run`, collect those prepared runs in a concretely
typed vector, and pass that vector here. Prepared runs accept no additional
`PASSVALUE` or execution keywords because their native Julia call shape is
already resolved.
"""
function prop_run_multi!(
    stack::AbstractArray{<:Number,3},
    samplings::AbstractVector{<:AbstractFloat},
    routine_name,
    lambda0_microns::Real,
    gridsize::Integer;
    PASSVALUE=nothing,
    kwargs...,
)
    prepared = prepare_prescription(
        routine_name,
        lambda0_microns,
        gridsize;
        PASSVALUE=PASSVALUE,
        kwargs...,
    )
    return prop_run_multi!(stack, samplings, prepared)
end

function prop_run_multi!(
    stack::AbstractArray{<:Number,3},
    samplings::AbstractVector{<:AbstractFloat},
    prepared::PreparedPrescription;
    PASSVALUE=prepared.passvalue,
    kwargs...,
)
    passes = _multi_passes(PASSVALUE)
    _validate_multi_run_destinations(stack, samplings, length(passes))
    batch = prepare_prescription_batch(prepared; pool_size=length(passes))
    return prop_run_multi!(stack, samplings, batch; PASSVALUE=passes, kwargs...)
end

function prop_run_multi!(
    stack::AbstractArray{<:Number,3},
    samplings::AbstractVector{<:AbstractFloat},
    batch::PreparedBatch;
    PASSVALUE=batch.prepared.passvalue,
    kwargs...,
)
    passes = _multi_passes(PASSVALUE)
    n = _validate_multi_run_destinations(stack, samplings, length(passes))
    contexts = ensure_prepared_batch_contexts!(batch, n)
    items = @view contexts[1:n]
    style = _multi_run_mutating_exec_style(_multi_run_exec_style(items), stack, samplings)
    runner = MultiRunPreparedRunner(batch.prepared, (; kwargs...))
    return _prop_run_multi_items!(style, stack, samplings, items, passes, runner)
end

function prop_run_multi!(
    stack::AbstractArray{<:Number,3},
    samplings::AbstractVector{<:AbstractFloat},
    model::PreparedModel;
    PASSVALUE=model.prepared.passvalue,
    kwargs...,
)
    passes = _multi_passes(PASSVALUE)
    n = _validate_multi_run_destinations(stack, samplings, length(passes))
    contexts = ensure_prepared_batch_contexts!(model.batch, n)
    model.assets isa PreparedAssetPool && ensure_prepared_asset_pool!(model.assets, n)
    style = _multi_run_mutating_exec_style(
        _multi_run_exec_style(@view contexts[1:n]),
        stack,
        samplings,
    )
    runner = MultiRunModelRunner(model, (; kwargs...))
    return _prop_run_multi_items!(style, stack, samplings, Base.OneTo(n), passes, runner)
end

function prop_run_multi!(
    stack::AbstractArray{<:Number,3},
    samplings::AbstractVector{<:AbstractFloat},
    runs::AbstractVector{<:MultiRunItem};
    PASSVALUE=nothing,
    kwargs...,
)
    passes = _multi_run_vector_passes(PASSVALUE, length(runs))
    _validate_multi_run_destinations(stack, samplings, length(runs))
    style = _multi_run_mutating_exec_style(_multi_run_exec_style(runs), stack, samplings)
    return _prop_run_multi_items!(
        style,
        stack,
        samplings,
        runs,
        passes,
        MultiRunItemRunner((; kwargs...)),
    )
end
