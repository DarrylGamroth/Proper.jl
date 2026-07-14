using BenchmarkTools
using Proper

const PREPARED_STEADY_GRID_N = 512
const PREPARED_STEADY_SAMPLES = 20
const PREPARED_BATCH_WAVELENGTHS = (0.50, 0.55, 0.60, 0.65)

function prepared_steady_state_prescription(λm, n; kwargs...)
    T = typeof(λm)
    wf = prop_begin(T(2.4), λm, n; beam_diam_fraction=T(0.5))
    prop_circular_aperture(wf, T(0.6))
    prop_lens(wf, T(20.0))
    prop_propagate(wf, T(20.0))
    return wf
end

@inline cpu_prepared_context(::Type{T}, grid_n::Integer) where {T<:AbstractFloat} = RunContext(Matrix{T})

function backend_prepared_context(beginfn::F, ::Type{T}, grid_n::Integer) where {F<:Function,T<:AbstractFloat}
    wf = beginfn(T, 2.4, 0.55e-6, grid_n; beam_diam_fraction=T(0.5))
    return RunContext(wf)
end

function prepare_steady_state_model(
    ::Type{T},
    grid_n::Integer,
    context::Union{Nothing,RunContext};
    name=:steady_state,
    wavelength_microns::Real=T(0.55),
) where {T<:AbstractFloat}
    return prepare_model(
        name,
        prepared_steady_state_prescription,
        wavelength_microns,
        grid_n;
        context=context,
        pool_size=1,
    )
end

function prepare_steady_state_sweep(
    ::Type{T},
    context_factory::F,
    wavelengths,
    grid_n::Integer=PREPARED_STEADY_GRID_N;
    name_prefix::Symbol=:steady_state_batch,
) where {T<:AbstractFloat,F<:Function}
    return [
        prepare_steady_state_model(
            T,
            grid_n,
            context_factory(T, grid_n);
            name=Symbol(name_prefix, :_, i),
            wavelength_microns=T(λm),
        ) for (i, λm) in enumerate(wavelengths)
    ]
end

function preallocated_steady_state_prescription(
    λm,
    n;
    wavefront,
    output,
    run_context,
)
    T = typeof(λm)
    prop_begin!(wavefront, T(2.4), λm; beam_diam_fraction=T(0.5))
    prop_circular_aperture(wavefront, T(0.6))
    prop_lens(wavefront, T(20), run_context)
    prop_propagate(wavefront, T(20), run_context)
    return prop_end(wavefront, output)
end

function prepare_preallocated_cpu_sweep(
    wavelengths,
    grid_n::Integer,
)
    T = Float64
    stack = zeros(T, grid_n, grid_n, length(wavelengths))
    samplings = zeros(T, length(wavelengths))
    runs = map(enumerate(wavelengths)) do (slot, wavelength_microns)
        context = cpu_prepared_context(T, grid_n)
        field = Matrix{Complex{T}}(undef, grid_n, grid_n)
        wavefront = prop_begin!(
            field,
            T(2.4),
            T(wavelength_microns) * T(1e-6);
            beam_diam_fraction=T(0.5),
            context=context,
        )
        prepared = prepare_prescription(
            preallocated_steady_state_prescription,
            T(wavelength_microns),
            grid_n;
            context=context,
        )
        return prepare_run(
            prepared;
            activate_context=false,
            wavefront=wavefront,
            output=@view(stack[:, :, slot]),
            run_context=context,
        )
    end
    return stack, samplings, runs
end

function run_prepared_workload(prepared, sync!::F=()->nothing) where {F<:Function}
    prop_run(prepared)
    sync!()
    return nothing
end

function run_prepared_sweep_workload(prepared_runs, sync!::F=()->nothing) where {F<:Function}
    prop_run_multi(prepared_runs)
    sync!()
    return nothing
end

function benchmark_prepared_trial(workload::F; samples::Integer=PREPARED_STEADY_SAMPLES) where {F<:Function}
    workload()
    workload()
    return run(@benchmarkable $workload() evals=1 samples=samples)
end
