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

@inline prepared_context(prepared::PreparedPrescription) = prepared.context

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
