using Random

struct RunContext{RNGT,B<:BackendStyle,FS<:FFTStyle,IS<:InterpStyle}
    rng::RNGT
    backend::B
    fft::FS
    interp::IS
    verbose::Bool
end

function RunContext(; rng=Random.default_rng(), verbose::Bool=false)
    return RunContext(Matrix; rng=rng, verbose=verbose)
end

function RunContext(::Type{A}; rng=Random.default_rng(), verbose::Bool=false) where {A<:AbstractArray}
    b = backend_style(A)
    f = fft_style(A)
    i = interp_style(A)
    return RunContext{typeof(rng),typeof(b),typeof(f),typeof(i)}(rng, b, f, i, verbose)
end
