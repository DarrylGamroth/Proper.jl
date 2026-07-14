# Prepared Execution Guide

## Purpose
This guide explains the prepared execution layer in `Proper.jl`:

- `PreparedPrescription`
- `PreparedBatch`
- `PreparedAssetPool`
- `PreparedModel`

These APIs are the recommended path when a prescription will be run more than
once and you want explicit ownership of runtime state, context reuse, or cached
assets.

## When To Use It
Use the plain compatibility surface first when you are porting or validating a
prescription:

```julia
psf, sampling = prop_run(my_prescription, 0.55, 256; dm=zeros(256, 256))
```

Move to the prepared layer when you need one or more of:

- repeated runs with the same prescription shape
- explicit `RunContext` reuse
- parallel repeated execution
- cached static assets
- a named execution object to hand around in application code

## Layering

### Upstream PROPER Translation
For users coming from upstream PROPER, the prepared execution layer should be
read as a progression from the familiar single-call surface:

- plain upstream-style execution:
  - `prop_run(my_prescription, 0.55, 256; use_dm=true)`
- same prescription reused:
  - `prepare_prescription(...)`
- same prescription reused in repeated/parallel passes:
  - `prepare_prescription_batch(...)`
- same prescription reused as one named application object with cached assets:
  - `prepare_model(...)`
- wavelength sweep:
  - build one prepared run per wavelength and call `prop_run_multi(runs)`

That is the intended Julia model. `PreparedModel` is the reusable execution
object; it is not a separate optical model type.

### `PreparedPrescription`
Use this when you want to normalize and retain a single prescription call shape.

```julia
prepared = prepare_prescription(my_prescription, 0.55, 256; use_dm=false)
psf, sampling = prop_run(prepared)
```

What it stores:

- resolved routine
- wavelength in meters
- grid size
- optional prepared `RunContext`
- optional explicit execution precision
- fixed kwargs
- default compatibility `PASSVALUE`

What it does not add:

- per-slot parallel contexts
- cached assets

Precision can be requested explicitly when you want the prepared execution
object to own a `Float32` or `Float64` runtime surface:

```julia
prepared = prepare_prescription(my_prescription, 0.55, 256; precision=Float32)
psf, sampling = prop_run(prepared)
```

If you also pass `context=...`, the context workspace precision must match the
requested `precision`.

### `PreparedBatch`
Use this when you want repeated or parallel execution with reusable per-slot
contexts.

```julia
prepared = prepare_prescription(my_prescription, 0.55, 256)
batch = prepare_prescription_batch(prepared; pool_size=4)

psf1, sampling1 = prop_run(batch; slot=1, use_dm=false)
psf2, sampling2 = prop_run(batch; slot=2, use_dm=true, dm=zeros(256, 256))
```

For parallel runs:

```julia
stack, samplings = prop_run_multi(
    batch;
    PASSVALUE=[
        (; use_dm=false),
        (; use_dm=true, dm=zeros(256, 256)),
    ],
)
```

`prop_run_multi` uses `PASSVALUE` for per-run variation. Map-like entries are
still normalized into native keywords before each prescription call.

Key behavior:

- slots are 1-based
- each slot gets an independent `RunContext`
- output order matches input order
- repeated runs preserve the prepared backend and precision

Reset reusable contexts with:

```julia
reset_prepared_batch!(batch)
```

This clears reusable workspace caches but intentionally does not rewind RNG
streams. Recreate the prepared contexts (or supply freshly seeded contexts) for
an exact stochastic replay.

## Coherent Carrier Phase

PROPER normally propagates the slowly varying complex envelope. For coherent
instrument arms such as an SCC reference beam, opt in to the uniform carrier
phase so path-length differences survive recombination:

```julia
prepared = prepare_prescription(
    my_prescription,
    0.55,
    512;
    phase_offset=true,
)
```

The equivalent typed context is:

```julia
ctx = RunContext(Matrix{Float64}; carrier_phase=TrackCarrierPhase())
```

`PHASE_OFFSET=true` is also accepted for upstream compatibility. The default is
`EnvelopeOnly()`, matching Python/MATLAB PROPER. Intensity from a single arm is
unchanged either way, so coherent simulations should validate the complex field
or the recombined-arm intensity.

## CPU Threading

PROPER's CPU propagation path is usually FFT-bound. FFTW and Julia task
threading are independent, and combining large counts can oversubscribe the
machine.

For one large serial prescription, a useful starting point is one Julia thread
and a modest FFTW count:

```julia
using Proper

prop_fftw_threads(4)  # call before constructing or warming prepared contexts
prepared = prepare_prescription(my_prescription, 0.55, 512)
psf, sampling = prop_run(prepared)
```

The equivalent process-level configuration is:

```sh
PROPER_FFTW_THREADS=4 julia --threads=1 --project=. my_run.jl
```

For wavelength sweeps and other `prop_run_multi` workloads, prefer outer Julia
parallelism with single-thread FFTW as the first configuration to benchmark:

```sh
PROPER_FFTW_THREADS=1 julia --threads=4 --project=. my_batch_run.jl
```

This outer threading applies to independent CPU contexts. AMDGPU, CUDA, and
unknown backends use ordered serial host submission by default so Julia tasks
do not contend while submitting whole prescriptions to a single device.
When a CPU batch already has enough independent runs to occupy the thread
pool, PROPER raises the crossover for inner cubic-interpolation and rotation
task fan-out. Smaller batches and standalone large interpolation or rotation
calls can use otherwise-idle Julia threads at smaller grid sizes.

Thread counts are machine- and grid-size-dependent. Benchmark representative
256, 512, and 1024 grids before choosing a production value. Changing
`prop_fftw_threads` affects newly created FFTW plans only; call
`reset_prepared_batch!` or recreate prepared contexts after changing it.

The repository includes a fresh-process matrix benchmark that rejects a
configuration if its complex field differs from the matrix reference before
timing it:

```sh
PROPER_THREAD_TOPOLOGY_WORKLOAD=core \
  PROPER_THREAD_TOPOLOGY_MATRIX="1:1 1:4 1:8 1:16" \
  PROPER_BENCH_GRID_N=512 scripts/benchmark_thread_topology.sh
```

Entries are `JULIA_THREADS:FFTW_THREADS`. On a Threadripper, include serial
outer/multithreaded FFTW points in the `core` workload. Then repeat with
`PROPER_THREAD_TOPOLOGY_WORKLOAD=batch`, including multithreaded
outer/single-thread FFTW and a few hybrid points whose product is near the
number of physical cores. Set `PROPER_BENCH_BATCH_SIZE` to the representative
number of independent runs; the default is four. Do not assume that all logical
threads are optimal: FFT scaling and memory bandwidth often flatten first,
while nested Julia and FFTW pools can oversubscribe the CPU. PROPER does not
choose or mutate a global thread topology automatically because the optimum
depends on grid size, prescription shape, concurrent runs, and other work in
the process.

The main propagation path is not BLAS-heavy. When Julia threads run independent
prescriptions, `LinearAlgebra.BLAS.set_num_threads(1)` is a sensible starting
point to avoid nested BLAS parallelism, but it should also be measured on the
application's non-propagation workloads.

### `PreparedAssetPool`
Use this when each slot should lazily create and retain its own asset bundle.

```julia
AssetBundle = @NamedTuple{dm::Matrix{Float64}, label::String}
assets = prepare_asset_pool(AssetBundle; pool_size=2) do slot
    return (dm=zeros(256, 256), label="slot_$slot")
end
```

The concrete asset type is required because a lazy factory may depend on the
eventual model. Declaring it keeps cache access type-stable without executing
the factory during setup. Factory results must be instances of that type;
mismatches raise `ArgumentError` instead of entering the cache.

Factory call forms accepted by the pool:

- `factory(slot, model)`
- `factory(slot)`
- `factory(model)`
- `factory()`

Factories for distinct slots may run concurrently during threaded batch
execution. Mutable state captured by a factory remains caller-owned and must be
thread-safe; the pool synchronizes neither that external state nor concurrent
reuse of the same slot.

Reset cached assets with:

```julia
reset_prepared_assets!(assets)
```

### `PreparedModel`
Use this when you want a named reusable execution object with both context reuse
and optional assets.

```julia
AssetBundle = @NamedTuple{dm::Matrix{Float64}, slot::Int}
model = prepare_model(
    :multi_example,
    multi_example,
    0.55,
    256;
    pool_size=2,
    assets=prepare_asset_pool(AssetBundle; pool_size=2) do slot
        return (dm=zeros(256, 256), slot=slot)
    end,
)
```

Single-slot execution:

```julia
psf, sampling = prop_run(model; slot=1, use_dm=true)
```

Parallel execution:

```julia
stack, samplings = prop_run_multi(
    model;
    PASSVALUE=[
        (; use_dm=false),
        (; use_dm=true),
    ],
)
```

Reset everything owned by the model with:

```julia
reset_prepared_model!(model)
```

You can also construct a vector of independently prepared runs and execute them
through the same stacked return surface:

```julia
runs = [
    prepare_model(:λ50, my_prescription, 0.50, 256; precision=Float32, pool_size=1),
    prepare_model(:λ55, my_prescription, 0.55, 256; precision=Float32, pool_size=1),
    prepare_model(:λ60, my_prescription, 0.60, 256; precision=Float32, pool_size=1),
]

stack, samplings = prop_run_multi(runs)
```

This is the intended prepared surface for wavelength sweeps or mixed prepared
execution objects. When the first output is already on a GPU backend,
`prop_run_multi` preserves that backend for the stacked output where feasible.

This vector form is also the throughput-oriented execution surface used by the
batch benchmark lane. It is the right API when your real workload is:

- fixed-shape wavelength sweeps
- repeated prepared runs
- GPU throughput rather than one-off interactive execution

Representative validated CUDA batch result for that surface:

- 4 prepared wavelengths on a `512 x 512` grid
- FP64: `7.82 ms`
- FP32: `537.67 us`

That is why explicit prepared precision selection and vector-of-prepared-runs
execution are documented as public capabilities rather than benchmark internals.

## Asset Injection Contract
When a `PreparedModel` resolves assets for a slot:

- if the asset value is a `NamedTuple`, its fields are merged directly into the
  execution kwargs
- otherwise it is passed as `assets=...`

That means both of these are valid:

```julia
AssetBundle = @NamedTuple{dm::Matrix{Float64}, use_dm::Bool}
assets = prepare_asset_pool(() -> (dm=zeros(256, 256), use_dm=true), AssetBundle; pool_size=1)
model = prepare_model(my_prescription, 0.55, 256; assets=assets, pool_size=1)
```

```julia
assets = prepare_asset_pool(
    () -> Dict("dm" => zeros(256, 256)),
    Dict{String,Matrix{Float64}};
    pool_size=1,
)
model = prepare_model(my_prescription, 0.55, 256; assets=assets, pool_size=1)
```

In the first case, `dm` and `use_dm` become direct kwargs. In the second case,
the prescription receives `assets=Dict(...)`.

## Recommended Usage Pattern

### Porting And Parity
Start here:

1. port the prescription with the familiar `prop_*` calls
2. run it with `prop_run(...)`
3. verify parity against the Python baseline

### Repeated Local Runs
Then move to:

1. `prepare_prescription(...)`
2. `prepare_prescription_batch(...)` if you need explicit slot reuse

### Application Or Model Layer
Use:

1. `prepare_model(...)`
2. `prepare_asset_pool(...)` if static assets should be cached per slot

## Relationship To The Compatibility Surface
The prepared APIs are not a second compatibility mode.

They sit under the familiar `prop_run` / `prop_run_multi` execution contract:

- plain calls are still supported
- prepared calls preserve the same return shapes
- prepared objects just make ownership, reuse, and batching explicit
- explicit `precision=Float32` / `Float64` is a numerics choice, not a separate
  execution mode

## Minimal Examples

### Single Prepared Prescription
```julia
prepared = prepare_prescription(my_prescription, 0.55, 256)
psf, sampling = prop_run(prepared)
```

### Batch With Parallel Passes
```julia
prepared = prepare_prescription(my_prescription, 0.55, 256)
batch = prepare_prescription_batch(prepared; pool_size=2)

stack, samplings = prop_run_multi(
    batch;
    PASSVALUE=[(; a=1), (; a=2)],
)
```

### Model With Cached Assets
```julia
AssetBundle = @NamedTuple{label::String, dm::Matrix{Float64}}
assets = prepare_asset_pool(AssetBundle; pool_size=2) do slot
    return (label="slot_$slot", dm=zeros(256, 256))
end

model = prepare_model(
    :my_model,
    my_prescription,
    0.55,
    256;
    assets=assets,
    pool_size=2,
)

stack, samplings = prop_run_multi(
    model;
    PASSVALUE=[(; use_dm=false), (; use_dm=true)],
)
```

### Precision-Selected Batch Sweep
```julia
runs = [
    prepare_prescription(my_prescription, 0.50, 256; precision=Float32),
    prepare_prescription(my_prescription, 0.55, 256; precision=Float32),
    prepare_prescription(my_prescription, 0.60, 256; precision=Float32),
]

stack, samplings = prop_run_multi(runs)
```

### Fully Prepared Single Runs
Use `prepare_run` only after the prescription call shape is stable and you need
to remove repeated model asset resolution, slot lookup, and keyword merging
from a repeated execution loop.

```julia
model = prepare_model(my_prescription, 0.55, 256; pool_size=1)
prepared = prepare_run(model; payload=payload)

psf, sampling = prop_run(prepared)
```

For most users, `prop_run(model; kwargs...)` remains the right prepared API.
`prop_run(prepared)` is the lower-level form for application loops where the
same mutable payload and output buffers are reused frame after frame.

`prepare_run` intentionally supports native Julia keyword prescriptions only.
`PASSVALUE` remains available through `prop_run` as the upstream compatibility
adapter, but it is not part of this fully prepared surface.

By default, prepared runs still activate the stored `RunContext` around the
prescription:

```julia
prepared = prepare_run(model; payload=payload)
```

Set `activate_context=false` only when the prescription explicitly passes the
stored context to context-aware operations. This avoids task-local context setup
in the hot loop, but the prescription must own that explicit context threading:

```julia
function hil_prescription(λm, n; payload, wavefront, output, run_context)
    wf = prop_begin!(wavefront, payload.diam_m, λm; beam_diam_fraction=1.0)
    prop_lens(wf, payload.focal_length_m, run_context, "science_lens")
    return prop_end(wf, output)
end

ctx = RunContext(Matrix{Float32})
model = prepare_model(hil_prescription, 0.55, 256; context=ctx, pool_size=1)
run_context = prepared_context(model, 1)
field = Matrix{ComplexF32}(undef, 256, 256)
wavefront = prop_begin!(field, payload.diam_m, 0.55f-6;
    beam_diam_fraction=1.0,
    context=run_context,
)
output = Matrix{Float32}(undef, 256, 256)
prepared = prepare_run(
    model;
    slot=1,
    activate_context=false,
    payload=payload,
    wavefront=wavefront,
    output=output,
    run_context=run_context,
)
```

This is the intended boundary for low-latency HIL and RTC integration: resolve
all static structure once, reuse caller-owned arrays, and keep per-frame work
limited to command updates plus propagation.

The same ownership model extends to independent batches with
`prop_run_multi!`. Bind each prepared run to its own wavefront, context, and
output slice, then reuse the stack and sampling vector:

```julia
function batch_hil_prescription(λm, n; wavefront, output, run_context)
    prop_begin!(wavefront, 2.4f0, λm; beam_diam_fraction=0.5f0)
    prop_lens(wavefront, 20f0, run_context)
    prop_propagate(wavefront, 20f0, run_context)
    return prop_end(wavefront, output)
end

wavelengths = Float32[0.50, 0.55, 0.60, 0.65]
n = 256
stack = zeros(Float32, n, n, length(wavelengths))
samplings = zeros(Float32, length(wavelengths))

runs = map(enumerate(wavelengths)) do (slot, wavelength)
    context = RunContext(Matrix{Float32})
    field = Matrix{ComplexF32}(undef, n, n)
    wavefront = prop_begin!(
        field,
        2.4f0,
        wavelength * 1f-6;
        beam_diam_fraction=0.5f0,
        context=context,
    )
    prepared = prepare_prescription(
        batch_hil_prescription,
        wavelength,
        n;
        context=context,
    )
    prepare_run(
        prepared;
        activate_context=false,
        wavefront=wavefront,
        output=@view(stack[:, :, slot]),
        run_context=context,
    )
end

prop_run_multi!(stack, samplings, runs)
```

The vector must remain concretely typed, and every run must own independent
mutable workspace. Dense strided CPU slices can execute on Julia threads.
Packed arrays such as `BitArray` and custom storage use serial execution because
distinct logical slices can share physical words. Accelerator batches retain
ordered serial host submission while their kernels execute on the device.

This is preferable to an arena or bump allocator for the current workload:
large arrays have stable shapes and lifetimes, so explicit typed buffers make
ownership, aliasing, backend placement, and numerical validation visible.
Arenas remain a poor fit for arrays that escape a prescription or must be
retained as batch output.
