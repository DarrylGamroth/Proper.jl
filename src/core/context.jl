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
    policy = resolve_compat_policy(compat_mode)
    b = backend_style(Array{Float64,2})
    f = fft_style(Array{Float64,2})
    i = interp_style(Array{Float64,2})
    return RunContext(policy, rng, b, f, i, verbose)
end

function RunContext(policy::P; rng=Random.default_rng(), verbose::Bool=false) where {P<:CompatPolicy}
    b = backend_style(Array{Float64,2})
    f = fft_style(Array{Float64,2})
    i = interp_style(Array{Float64,2})
    return RunContext{P,typeof(rng),typeof(b),typeof(f),typeof(i)}(policy, rng, b, f, i, verbose)
end
