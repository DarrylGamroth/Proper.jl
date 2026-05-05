# Adaptive Optics Integration Guide

Status: active

## Purpose

This guide describes the recommended boundary between `Proper.jl` and
AO/runtime packages such as `AdaptiveOpticsSim.jl`.

The intended HIL use case is:

1. an RTC provides actuator commands
2. an AO runtime converts those commands into a sampled OPD or DM surface
3. `Proper.jl` evaluates the science prescription and returns pixels

`Proper.jl` should remain the wave-optics propagation package. It should not
become the primary actuator-reconstruction runtime when a dedicated AO package
already owns that model.

## Recommended Boundary

For HIL and coronagraph science-arm integrations, pass sampled arrays into a
prepared Proper model through ordinary Julia keywords:

```julia
struct CoronagraphPayload{O,P,T}
    opd_m::O
    pupil::P
    diameter_m::T
    focal_length_m::T
    lyot_stop_norm::T
end

function coronagraph_prescription(λm, n; payload::CoronagraphPayload)
    wf = prop_begin(payload.diameter_m, λm, n; beam_diam_fraction=1.0)
    prop_multiply(wf, payload.pupil)
    prop_add_phase(wf, payload.opd_m)
    return prop_end(wf)
end

model = prepare_model(:science_arm, coronagraph_prescription, 1.65, 256; pool_size=1)
psf, sampling = prop_run(model; payload=payload)
```

Use `PASSVALUE` only when preserving upstream PROPER calling conventions or an
existing parity harness. New Julia-native integrations should prefer typed
payloads and normal keywords.

## OPD Versus DM Surface

There are two valid handoff styles:

- **Total OPD handoff:** the AO package has already applied atmosphere,
  telescope, and controllable-optic state. Use `prop_add_phase(wf, opd_m)`.
- **DM-surface handoff:** the Proper prescription intentionally owns the DM
  surface application. Use `prop_dm(wf, dm_surface)` when `dm_surface` is
  already sampled on the wavefront grid.

Do not put the upstream-compatible actuator-space `prop_dm` adapter in the HIL
hot loop unless the Proper prescription must own actuator geometry and
influence-function modeling. For RTC workflows, actuator-to-surface conversion
usually belongs in the AO runtime or application-specific reconstructor.

## Validation Required Before Support Claims

Before treating an AO/Proper integration as supported, add deterministic checks
for:

- OPD units in meters
- sign convention using piston, tilt, focus, and one actuator poke
- centering convention between the AO grid and the Proper wavefront grid
- row/column orientation using asymmetric maps
- pupil ownership, including central obstruction, spiders, segment gaps, and
  reflectivity when relevant
- same-backend behavior on CUDA and AMDGPU when GPU support is claimed

These checks should exist before relying on throughput numbers. Performance is
not meaningful if the map is transposed, shifted, or signed incorrectly.

## GPU Guidance

For same-device execution:

```julia
ctx = RunContext(typeof(opd_m))
model = prepare_model(:science_arm, coronagraph_prescription, 1.65, 256;
    context=ctx,
    precision=Float32,
    pool_size=1,
)
psf, sampling = prop_run(model; payload=payload)
```

The important GPU paths for this boundary are:

- `prop_add_phase` for sampled OPD maps
- `prop_dm(wf, dm_surface)` for already-sampled DM maps
- pupil multiplication and coronagraph masks
- prepared propagation and `prop_end`
- batched or multi-wavelength prepared execution when throughput matters

FITS IO and other file-backed map loading remain host boundaries. Transfer to a
GPU array explicitly before entering the repeated frame loop.

## Benchmark Guidance

Benchmark the combined command-to-pixels path, not only isolated kernels:

1. stage/update the AO command
2. update the sampled OPD or DM surface
3. run the prepared Proper science prescription
4. return the detector or science pixels needed by the RTC

Use `BenchmarkTools` after warmup and record allocations separately from
Julia cold-start / TTFx. When CUDA or AMDGPU support is claimed, validate on
real hardware rather than inferring support from CPU tests.

## Relationship To `AdaptiveOpticsSim.jl`

The recommended `AdaptiveOpticsSim.jl` direction is:

- `AdaptiveOpticsSim.jl` owns AO runtime state, command staging, WFS products,
  and DM-to-OPD conversion.
- `Proper.jl` owns external science prescriptions that need PROPER-compatible
  propagation or coronagraph elements.
- The application boundary passes typed payloads containing OPD, pupil, and
  scalar instrument parameters.

Keep this as an integration seam unless repeated use proves that a separate
small integration package is justified.
