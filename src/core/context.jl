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

struct RunContext{
    RNGT,
    B<:BackendStyle,
    FS<:FFTStyle,
    IS<:InterpStyle,
    PS<:FFTPlanningStyle,
    CS<:CarrierPhaseStyle,
    WS,
}
    rng::RNGT
    backend::B
    fft::FS
    interp::IS
    planning::PS
    carrier_phase::CS
    workspace::WS
    verbose::Bool
end

function RunContext(
    ;
    rng=Random.default_rng(),
    verbose::Bool=false,
    fft_planning::FFTPlanningStyle=FFTEstimateStyle(),
    carrier_phase::CarrierPhaseStyle=EnvelopeOnly(),
)
    return RunContext(
        Matrix;
        rng=rng,
        verbose=verbose,
        fft_planning=fft_planning,
        carrier_phase=carrier_phase,
    )
end

function RunContext(
    ::Type{A},
    ws::WS;
    rng=Random.default_rng(),
    verbose::Bool=false,
    fft_planning::FFTPlanningStyle=FFTEstimateStyle(),
    carrier_phase::CarrierPhaseStyle=EnvelopeOnly(),
) where {A<:AbstractArray,WS}
    b = backend_style(A)
    f = fft_style(A)
    i = interp_style(A)
    return RunContext{
        typeof(rng),
        typeof(b),
        typeof(f),
        typeof(i),
        typeof(fft_planning),
        typeof(carrier_phase),
        typeof(ws),
    }(
        rng,
        b,
        f,
        i,
        fft_planning,
        carrier_phase,
        ws,
        verbose,
    )
end

function RunContext(
    ::Type{A};
    rng=Random.default_rng(),
    verbose::Bool=false,
    fft_planning::FFTPlanningStyle=FFTEstimateStyle(),
    carrier_phase::CarrierPhaseStyle=EnvelopeOnly(),
) where {A<:AbstractArray}
    T = _workspace_float_type(A)
    return RunContext(
        A,
        ProperWorkspace(A, T);
        rng=rng,
        verbose=verbose,
        fft_planning=fft_planning,
        carrier_phase=carrier_phase,
    )
end

@inline function RunContext(
    wf::WaveFront;
    rng=Random.default_rng(),
    verbose::Bool=false,
    fft_planning::FFTPlanningStyle=FFTEstimateStyle(),
    carrier_phase::CarrierPhaseStyle=EnvelopeOnly(),
)
    return RunContext(
        typeof(wf.field),
        wf.workspace;
        rng=rng,
        verbose=verbose,
        fft_planning=fft_planning,
        carrier_phase=carrier_phase,
    )
end

@inline active_run_context() =
    get(task_local_storage(), ACTIVE_RUN_CONTEXT_KEY, nothing)::Union{Nothing,RunContext}
@inline active_run_workspace() = active_run_workspace(active_run_context())
@inline active_run_workspace(::Nothing) = nothing
@inline active_run_workspace(ctx::RunContext) = ctx.workspace

@inline _resolve_run_context(wf::WaveFront, ::Nothing) = RunContext(wf)
@inline _resolve_run_context(wf::WaveFront, ctx::RunContext) =
    wf.workspace === ctx.workspace ? ctx : RunContext(wf)
@inline resolve_run_context(wf::WaveFront) = _resolve_run_context(wf, active_run_context())

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
@inline carrier_phase_style(ctx::RunContext) = ctx.carrier_phase
@inline interp_workspace(ctx::RunContext) = ctx.workspace.interp
@inline mask_workspace(ctx::RunContext) = ctx.workspace.mask
@inline fft_workspace(ctx::RunContext) = ctx.workspace.fft
@inline sampling_workspace(ctx::RunContext) = ctx.workspace.sampling
@inline qphase_workspace(ctx::RunContext) = ctx.workspace.qphase

@inline function fresh_context(ctx::RunContext; rng=fork_rng(ctx.rng, 1))
    ws = fresh_workspace(ctx.workspace)
    A = workspace_backend_array_type(ws)
    return RunContext(
        A,
        ws;
        rng=rng,
        verbose=ctx.verbose,
        fft_planning=ctx.planning,
        carrier_phase=ctx.carrier_phase,
    )
end

@inline fresh_context(::Nothing; rng=Random.default_rng()) = nothing

@inline function with_carrier_phase(ctx::RunContext, carrier_phase::CarrierPhaseStyle)
    A = workspace_backend_array_type(ctx.workspace)
    return RunContext(
        A,
        ctx.workspace;
        rng=ctx.rng,
        verbose=ctx.verbose,
        fft_planning=ctx.planning,
        carrier_phase=carrier_phase,
    )
end

@inline _carrier_phase_style(enabled::Bool) = enabled ? TrackCarrierPhase() : EnvelopeOnly()

@inline function context_with_phase_override(
    ctx::Union{Nothing,RunContext},
    ::Nothing,
    ::Type{T},
) where {T<:AbstractFloat}
    return ctx
end

@inline function context_with_phase_override(
    ::Nothing,
    enabled::Bool,
    ::Type{T},
) where {T<:AbstractFloat}
    return RunContext(Matrix{T}; carrier_phase=_carrier_phase_style(enabled))
end

@inline function context_with_phase_override(
    ctx::RunContext,
    enabled::Bool,
    ::Type{T},
) where {T<:AbstractFloat}
    return with_carrier_phase(ctx, _carrier_phase_style(enabled))
end

@inline _apply_carrier_phase!(wf::WaveFront, ::Real, ::EnvelopeOnly) = wf

@inline function _apply_carrier_phase!(
    wf::WaveFront{T},
    dz::Real,
    ::TrackCarrierPhase,
) where {T<:AbstractFloat}
    P = promote_type(T, typeof(float(dz)))
    wavelength = P(wf.wavelength_m)
    reduced_distance = rem(P(dz), wavelength, RoundNearest)
    factor = Complex{T}(cispi(P(2) * reduced_distance / wavelength))
    wf.field .*= factor
    return wf
end

@inline apply_carrier_phase!(wf::WaveFront, dz::Real, ctx::RunContext) =
    _apply_carrier_phase!(wf, dz, carrier_phase_style(ctx))
