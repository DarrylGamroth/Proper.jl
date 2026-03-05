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

function RunContext(::Type{A}; rng=Random.default_rng(), verbose::Bool=false) where {A<:AbstractArray}
    b = backend_style(A)
    f = fft_style(A)
    i = interp_style(A)
    T = _workspace_float_type(A)
    ws = ProperWorkspace(T)
    return RunContext{typeof(rng),typeof(b),typeof(f),typeof(i),typeof(ws)}(rng, b, f, i, ws, verbose)
end

@inline fft_style(ctx::RunContext) = ctx.fft
@inline interp_style(ctx::RunContext) = ctx.interp
@inline interp_workspace(ctx::RunContext) = ctx.workspace.interp
@inline fft_workspace(ctx::RunContext) = ctx.workspace.fft
