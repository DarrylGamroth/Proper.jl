"""Run a prescription function by function object or global name."""
function prop_run(
    routine_name,
    lambda0_microns::Real,
    gridsize::Integer;
    PASSVALUE=nothing,
    context::Union{Nothing,RunContext}=nothing,
    kwargs...,
)
    λm = float(lambda0_microns) * 1e-6
    fn = if routine_name isa Function
        routine_name
    elseif routine_name isa Symbol
        getfield(Main, routine_name)
    elseif routine_name isa AbstractString
        getfield(Main, Symbol(routine_name))
    else
        throw(ArgumentError("routine_name must be Function, Symbol, or String"))
    end

    result = with_run_context(context) do
        if PASSVALUE === nothing
            fn(λm, gridsize; kwargs...)
        else
            try
                fn(λm, gridsize, PASSVALUE; kwargs...)
            catch err
                if err isa MethodError
                    fn(λm, gridsize; PASSVALUE=PASSVALUE, kwargs...)
                else
                    rethrow(err)
                end
            end
        end
    end

    if result isa WaveFront
        return prop_end(result)
    elseif result isa Tuple && length(result) == 2
        return result
    end
    throw(ArgumentError("Prescription must return WaveFront or (psf, sampling)"))
end
