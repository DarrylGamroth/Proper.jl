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
psf, sampling = prop_run(my_prescription, 0.55, 256; PASSVALUE=Dict("dm" => zeros(48, 48)))
```

Move to the prepared layer when you need one or more of:

- repeated runs with the same prescription shape
- explicit `RunContext` reuse
- parallel repeated execution
- cached static assets
- a named execution object to hand around in application code

## Layering

### `PreparedPrescription`
Use this when you want to normalize and retain a single prescription call shape.

```julia
prepared = prepare_prescription(my_prescription, 0.55, 256; PASSVALUE=Dict("use_dm" => false))
psf, sampling = prop_run(prepared)
```

What it stores:

- resolved routine
- wavelength in meters
- grid size
- optional prepared `RunContext`
- fixed kwargs
- default `PASSVALUE`

What it does not add:

- per-slot parallel contexts
- cached assets

### `PreparedBatch`
Use this when you want repeated or parallel execution with reusable per-slot
contexts.

```julia
prepared = prepare_prescription(my_prescription, 0.55, 256)
batch = prepare_prescription_batch(prepared; pool_size=4)

psf1, sampling1 = prop_run(batch; slot=1, PASSVALUE=Dict("use_dm" => false))
psf2, sampling2 = prop_run(batch; slot=2, PASSVALUE=Dict("use_dm" => true, "dm" => zeros(48, 48)))
```

For parallel runs:

```julia
stack, samplings = prop_run_multi(
    batch;
    PASSVALUE=[
        Dict("use_dm" => false),
        Dict("use_dm" => true, "dm" => zeros(48, 48)),
    ],
)
```

Key behavior:

- slots are 1-based
- each slot gets an independent `RunContext`
- output order matches input order

Reset reusable contexts with:

```julia
reset_prepared_batch!(batch)
```

### `PreparedAssetPool`
Use this when each slot should lazily create and retain its own asset bundle.

```julia
assets = prepare_asset_pool(; pool_size=2) do slot
    return (dm=zeros(48, 48), label="slot_$slot")
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
        return (dm=zeros(48, 48), slot=slot)
    end,
)
```

Single-slot execution:

```julia
psf, sampling = prop_run(model; slot=1, PASSVALUE=Dict("use_dm" => true))
```

Parallel execution:

```julia
stack, samplings = prop_run_multi(
    model;
    PASSVALUE=[
        Dict("use_dm" => false),
        Dict("use_dm" => true),
    ],
)
```

Reset everything owned by the model with:

```julia
reset_prepared_model!(model)
```

## Asset Injection Contract
When a `PreparedModel` resolves assets for a slot:

- if the asset value is a `NamedTuple`, its fields are merged directly into the
  execution kwargs
- otherwise it is passed as `assets=...`

That means both of these are valid:

```julia
assets = prepare_asset_pool(() -> (dm=zeros(48, 48), use_dm=true); pool_size=1)
model = prepare_model(my_prescription, 0.55, 256; assets=assets, pool_size=1)
```

```julia
assets = prepare_asset_pool(() -> Dict("dm" => zeros(48, 48)); pool_size=1)
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
    PASSVALUE=[Dict("a" => 1), Dict("a" => 2)],
)
```

### Model With Cached Assets
```julia
assets = prepare_asset_pool(; pool_size=2) do slot
    return (label="slot_$slot", dm=zeros(48, 48))
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
    PASSVALUE=[Dict("use_dm" => false), Dict("use_dm" => true)],
)
```
