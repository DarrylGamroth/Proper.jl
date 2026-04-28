using Random

const ACTIVE_RUN_CONTEXT_KEY = :Proper_ACTIVE_RUN_CONTEXT

@inline function _seed_child_rng!(rng, token::Integer)
    return xor(rand(rng, UInt64), UInt64(token) * 0x9e3779b97f4a7c15)
end

@inline function fork_rng(rng, token::Integer)
    return Random.Xoshiro(_seed_child_rng!(rng, token))
end

@inline function _workspace_float_type(::Type{A}) where {A<:AbstractArray}
    ET = eltype(A)
    if ET <: Number
        return float(real(ET))
    end
    return Float64
end

struct RunContext{RNGT,B<:BackendStyle,FS<:FFTStyle,IS<:InterpStyle,PS<:FFTPlanningStyle,WS}
    rng::RNGT
    backend::B
    fft::FS
    interp::IS
    planning::PS
    workspace::WS
    verbose::Bool
end

function RunContext(; rng=Random.default_rng(), verbose::Bool=false, fft_planning::FFTPlanningStyle=FFTEstimateStyle())
    return RunContext(Matrix; rng=rng, verbose=verbose, fft_planning=fft_planning)
end

function RunContext(
    ::Type{A},
    ws::WS;
    rng=Random.default_rng(),
    verbose::Bool=false,
    fft_planning::FFTPlanningStyle=FFTEstimateStyle(),
) where {A<:AbstractArray,WS}
    b = backend_style(A)
    f = fft_style(A)
    i = interp_style(A)
    return RunContext{typeof(rng),typeof(b),typeof(f),typeof(i),typeof(fft_planning),typeof(ws)}(
        rng,
        b,
        f,
        i,
        fft_planning,
        ws,
        verbose,
    )
end

function RunContext(
    ::Type{A};
    rng=Random.default_rng(),
    verbose::Bool=false,
    fft_planning::FFTPlanningStyle=FFTEstimateStyle(),
) where {A<:AbstractArray}
    T = _workspace_float_type(A)
    return RunContext(A, ProperWorkspace(A, T); rng=rng, verbose=verbose, fft_planning=fft_planning)
end

@inline function RunContext(
    wf::WaveFront;
    rng=Random.default_rng(),
    verbose::Bool=false,
    fft_planning::FFTPlanningStyle=FFTEstimateStyle(),
)
    return RunContext(typeof(wf.field), wf.workspace; rng=rng, verbose=verbose, fft_planning=fft_planning)
end

@inline active_run_context() = get(task_local_storage(), ACTIVE_RUN_CONTEXT_KEY, nothing)
@inline active_run_workspace() = active_run_workspace(active_run_context())
@inline active_run_workspace(::Nothing) = nothing
@inline active_run_workspace(ctx::RunContext) = ctx.workspace

@inline resolve_run_workspace(::Nothing, ::Nothing) = active_run_workspace()
@inline resolve_run_workspace(ctx::RunContext, ::Nothing) = ctx.workspace
@inline resolve_run_workspace(::Nothing, ws::ProperWorkspace) = ws

@inline function resolve_run_workspace(ctx::RunContext, ws::ProperWorkspace)
    ws === ctx.workspace || throw(ArgumentError("workspace must match context.workspace"))
    return ws
end

@inline with_run_context(f::F, ::Nothing) where {F<:Function} = f()

@inline function with_run_context(f::F, ctx::RunContext) where {F<:Function}
    return task_local_storage(ACTIVE_RUN_CONTEXT_KEY, ctx) do
        f()
    end
end

@inline fft_style(ctx::RunContext) = ctx.fft
@inline interp_style(ctx::RunContext) = ctx.interp
@inline fft_planning_style(ctx::RunContext) = ctx.planning
@inline interp_workspace(ctx::RunContext) = ctx.workspace.interp
@inline mask_workspace(ctx::RunContext) = ctx.workspace.mask
@inline fft_workspace(ctx::RunContext) = ctx.workspace.fft
@inline sampling_workspace(ctx::RunContext) = ctx.workspace.sampling
@inline qphase_workspace(ctx::RunContext) = ctx.workspace.qphase

@inline function fresh_context(ctx::RunContext; rng=fork_rng(ctx.rng, 1))
    ws = fresh_workspace(ctx.workspace)
    A = workspace_backend_array_type(ws)
    return RunContext(A, ws; rng=rng, verbose=ctx.verbose, fft_planning=ctx.planning)
end

@inline fresh_context(::Nothing; rng=Random.default_rng()) = nothing
