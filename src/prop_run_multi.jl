using Base.Threads

abstract type MultiRunExecStyle end
struct MultiRunThreadedExecStyle <: MultiRunExecStyle end
struct MultiRunSerialExecStyle <: MultiRunExecStyle end

@inline _multi_run_exec_style(::CPUBackend) = MultiRunThreadedExecStyle()
@inline _multi_run_exec_style(::BackendStyle) = MultiRunSerialExecStyle()
@inline _multi_run_exec_style(::Nothing) = MultiRunThreadedExecStyle()
@inline _multi_run_exec_style(ctx::RunContext) = _multi_run_exec_style(ctx.backend)

@inline _multi_run_context(prepared::PreparedPrescription) = prepared.context
@inline _multi_run_context(batch::PreparedBatch) = prepared_context(batch, 1)
@inline _multi_run_context(model::PreparedModel) = prepared_context(model, 1)

function _multi_run_exec_style(contexts::AbstractVector)
    workspaces = IdDict{Any,Nothing}()
    for ctx in contexts
        _multi_run_exec_style(ctx) isa MultiRunSerialExecStyle && return MultiRunSerialExecStyle()
        if ctx isa RunContext
            ws = ctx.workspace
            haskey(workspaces, ws) && return MultiRunSerialExecStyle()
            workspaces[ws] = nothing
        end
    end
    return MultiRunThreadedExecStyle()
end

function _multi_run_exec_style(
    runs::AbstractVector{<:Union{PreparedPrescription,PreparedBatch,PreparedModel}},
)
    contexts = map(_multi_run_context, runs)
    return _multi_run_exec_style(contexts)
end

@inline _multi_passes(passvalue) = passvalue === nothing ? [nothing] : (passvalue isa AbstractVector ? passvalue : [passvalue])

@inline function _allocate_output_stack(first_out::AbstractMatrix, sy::Integer, sx::Integer, n::Integer)
    return similar(first_out, eltype(first_out), sy, sx, n)
end

@inline function _store_output_slice!(stack, i::Integer, out::AbstractMatrix)
    copyto!(selectdim(stack, 3, i), out)
    return stack
end

function _prop_run_multi_items(
    ::MultiRunSerialExecStyle,
    items::AbstractVector,
    passes,
    runner::F,
) where {F<:Function}
    n = length(items)
    n == length(passes) || throw(ArgumentError("PASSVALUE length must match number of prepared runs"))
    n > 0 || throw(ArgumentError("prepared run collection must not be empty"))

    first_out, first_sampling = runner(items[1], passes[1], 1)
    first_out isa AbstractMatrix || throw(ArgumentError("prop_run_multi expects matrix outputs"))

    outT = typeof(first_out)
    sampT = typeof(float(first_sampling))
    sy, sx = size(first_out)
    stack = _allocate_output_stack(first_out, sy, sx, n)
    samplings = Vector{sampT}(undef, n)
    _store_output_slice!(stack, 1, first_out)
    samplings[1] = float(first_sampling)

    @inbounds for i in 2:n
        out, s = runner(items[i], passes[i], i)
        out isa outT || throw(ArgumentError("All prop_run_multi outputs must have the same type"))
        size(out) == (sy, sx) || throw(ArgumentError("All prop_run_multi outputs must have the same size"))
        _store_output_slice!(stack, i, out)
        samplings[i] = float(s)
    end
    return stack, samplings
end

function _prop_run_multi_items(
    ::MultiRunThreadedExecStyle,
    items::AbstractVector,
    passes,
    runner::F,
) where {F<:Function}
    n = length(items)
    n == length(passes) || throw(ArgumentError("PASSVALUE length must match number of prepared runs"))
    n > 0 || throw(ArgumentError("prepared run collection must not be empty"))

    if Threads.nthreads() == 1 || n == 1
        return _prop_run_multi_items(MultiRunSerialExecStyle(), items, passes, runner)
    end

    results = Vector{Any}(undef, n)
    inner_min_elems = n >= Threads.nthreads() ? 262_144 : 0
    inner_policy = CPUInnerKernelPolicy(true, inner_min_elems)
    @threads for i in eachindex(items, passes)
        results[i] = with_cpu_inner_kernel_policy(inner_policy) do
            runner(items[i], passes[i], i)
        end
    end

    first_out, first_sampling = results[1]
    first_out isa AbstractMatrix || throw(ArgumentError("prop_run_multi expects matrix outputs"))

    outT = typeof(first_out)
    sampT = typeof(float(first_sampling))
    sy, sx = size(first_out)
    stack = _allocate_output_stack(first_out, sy, sx, n)
    samplings = Vector{sampT}(undef, n)

    @inbounds for i in eachindex(results)
        out, s = results[i]
        out isa outT || throw(ArgumentError("All prop_run_multi outputs must have the same type"))
        size(out) == (sy, sx) || throw(ArgumentError("All prop_run_multi outputs must have the same size"))
        _store_output_slice!(stack, i, out)
        samplings[i] = float(s)
    end

    return stack, samplings
end

"""
    prop_run_multi(routine_name, lambda0_microns, gridsize; PASSVALUE=nothing, kwargs...)
    prop_run_multi(prepared::PreparedPrescription; PASSVALUE=prepared.passvalue, kwargs...)
    prop_run_multi(batch::PreparedBatch; PASSVALUE=batch.prepared.passvalue, kwargs...)
    prop_run_multi(model::PreparedModel; PASSVALUE=model.prepared.passvalue, kwargs...)
    prop_run_multi(runs::AbstractVector{<:Union{PreparedPrescription,PreparedBatch,PreparedModel}}; PASSVALUE=nothing, kwargs...)

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
    style = _multi_run_exec_style(@view contexts[1:n])
    return _prop_run_multi_items(style, Base.OneTo(n), passes, (_, pass, slot) -> begin
        prop_run(model; PASSVALUE=pass, slot=slot, kwargs...)
    end)
end

function prop_run_multi(
    runs::AbstractVector{<:Union{PreparedPrescription,PreparedBatch,PreparedModel}};
    PASSVALUE=nothing,
    kwargs...,
)
    passes = PASSVALUE === nothing ? fill(nothing, length(runs)) :
        (PASSVALUE isa AbstractVector ? PASSVALUE : fill(PASSVALUE, length(runs)))
    style = _multi_run_exec_style(runs)
    return _prop_run_multi_items(style, runs, passes, (run, pass, _) -> begin
        prop_run(run; PASSVALUE=pass, kwargs...)
    end)
end
