using Random

struct RunContext{P<:CompatPolicy,RNGT,B<:BackendStyle,FS<:FFTStyle,IS<:InterpStyle}
    policy::P
    rng::RNGT
    backend::B
    fft::FS
    interp::IS
    verbose::Bool
end

function RunContext(; compat_mode::Symbol=:python334, rng=Random.default_rng(), verbose::Bool=false)
    return RunContext(Array{Float64,2}; compat_mode=compat_mode, rng=rng, verbose=verbose)
end

function RunContext(policy::P; rng=Random.default_rng(), verbose::Bool=false) where {P<:CompatPolicy}
    return RunContext(policy, Array{Float64,2}; rng=rng, verbose=verbose)
end

function RunContext(::Type{A}; compat_mode::Symbol=:python334, rng=Random.default_rng(), verbose::Bool=false) where {A<:AbstractArray}
    policy = resolve_compat_policy(compat_mode)
    return RunContext(policy, A; rng=rng, verbose=verbose)
end

function RunContext(policy::P, ::Type{A}; rng=Random.default_rng(), verbose::Bool=false) where {P<:CompatPolicy,A<:AbstractArray}
    b = backend_style(A)
    f = fft_style(A)
    i = interp_style(A)
    return RunContext{P,typeof(rng),typeof(b),typeof(f),typeof(i)}(policy, rng, b, f, i, verbose)
end
