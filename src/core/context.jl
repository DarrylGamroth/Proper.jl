using Random

@inline function _workspace_float_type(::Type{A}) where {A<:AbstractArray}
    ET = eltype(A)
    if ET <: Number
        return float(real(ET))
    end
    return Float64
end

struct RunContext{RNGT,B<:BackendStyle,FS<:FFTStyle,IS<:InterpStyle,WS}
    rng::RNGT
    backend::B
    fft::FS
    interp::IS
    workspace::WS
    verbose::Bool
end

function RunContext(; rng=Random.default_rng(), verbose::Bool=false)
    return RunContext(Matrix; rng=rng, verbose=verbose)
end

function RunContext(::Type{A}, ws::WS; rng=Random.default_rng(), verbose::Bool=false) where {A<:AbstractArray,WS}
    b = backend_style(A)
    f = fft_style(A)
    i = interp_style(A)
    return RunContext{typeof(rng),typeof(b),typeof(f),typeof(i),typeof(ws)}(rng, b, f, i, ws, verbose)
end

function RunContext(::Type{A}; rng=Random.default_rng(), verbose::Bool=false) where {A<:AbstractArray}
    T = _workspace_float_type(A)
    return RunContext(A, ProperWorkspace(A, T); rng=rng, verbose=verbose)
end

@inline function RunContext(wf::WaveFront; rng=Random.default_rng(), verbose::Bool=false)
    return RunContext(typeof(wf.field), wf.workspace; rng=rng, verbose=verbose)
end

@inline fft_style(ctx::RunContext) = ctx.fft
@inline interp_style(ctx::RunContext) = ctx.interp
@inline interp_workspace(ctx::RunContext) = ctx.workspace.interp
@inline mask_workspace(ctx::RunContext) = ctx.workspace.mask
@inline fft_workspace(ctx::RunContext) = ctx.workspace.fft
