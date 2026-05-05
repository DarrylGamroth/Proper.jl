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

### `PreparedAssetPool`
Use this when each slot should lazily create and retain its own asset bundle.

```julia
assets = prepare_asset_pool(; pool_size=2) do slot
    return (dm=zeros(256, 256), label="slot_$slot")
end
```

Factory call forms accepted by the pool:

- `factory(slot, model)`
- `factory(slot)`
- `factory(model)`
- `factory()`

Reset cached assets with:

```julia
reset_prepared_assets!(assets)
```

### `PreparedModel`
Use this when you want a named reusable execution object with both context reuse
and optional assets.

```julia
model = prepare_model(
    :multi_example,
    multi_example,
    0.55,
    256;
    pool_size=2,
    assets=prepare_asset_pool(; pool_size=2) do slot
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
assets = prepare_asset_pool(() -> (dm=zeros(256, 256), use_dm=true); pool_size=1)
model = prepare_model(my_prescription, 0.55, 256; assets=assets, pool_size=1)
```

```julia
assets = prepare_asset_pool(() -> Dict("dm" => zeros(256, 256)); pool_size=1)
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
assets = prepare_asset_pool(; pool_size=2) do slot
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

### Lower-Level Hot Calls
Use `prepare_hot_call` only after the prescription call shape is stable and you
need to remove repeated model asset resolution, slot lookup, and keyword
merging from a hot loop.

```julia
model = prepare_model(my_prescription, 0.55, 256; pool_size=1)
hot = prepare_hot_call(model; payload=payload)

psf, sampling = prop_run_hot(hot)
```

For most users, `prop_run(model; kwargs...)` remains the right prepared API.
`prop_run_hot(hot)` is a lower-level execution primitive for application loops
where the same mutable payload and output buffers are reused frame after frame.

`prepare_hot_call` intentionally supports native Julia keyword prescriptions
only. `PASSVALUE` remains available through `prop_run` as the upstream
compatibility adapter, but it is not part of this lower-level hot-loop surface.

By default, hot calls still activate the stored `RunContext` around the
prescription:

```julia
hot = prepare_hot_call(model; payload=payload)
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
hot = prepare_hot_call(
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
