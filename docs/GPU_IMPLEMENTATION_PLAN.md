# GPU Implementation Status

This document is the current maintainer tracker for GPU-facing work in
`Proper.jl`. Historical slice notes are archived in
[`archive/GPU_IMPLEMENTATION_PLAN_2026-03-17.md`](archive/GPU_IMPLEMENTATION_PLAN_2026-03-17.md).

## Current Contract
- CUDA and AMDGPU are optional package extensions.
- GPU performance work targets shared `Proper.jl` propagation and helper paths,
  not WFIRST-specific reference-model code.
- The intended performance surface is explicit `RunContext` plus mutating `!`
  APIs after warmup.
- Allocating convenience wrappers remain valid user APIs, but are not the
  zero-allocation benchmark target.
- FITS I/O remains host-backed; callers should treat host/device transfer at I/O
  boundaries as explicit.

## Completed Work
- backend traits route CPU, CUDA, AMDGPU, FFT, interpolation, and shift behavior
- CUDA and AMDGPU extensions provide backend-native propagation workspaces
- warmed propagation allocation gates cover core GPU propagation operations
- hidden host fallback was removed from GPU-visible hot paths that previously
  materialized through `Matrix(...)`
- CPU/CUDA/AMDGPU benchmark lanes use warmed steady-state `BenchmarkTools`
  measurements and report skipped GPU backends instead of failing when hardware
  is absent
- prepared vector-of-runs execution supports multi-wavelength throughput
- explicit `Float32` prepared execution is documented and benchmarked

## Active Follow-Up Candidates
- tighten GPU allocation thresholds as CUDA.jl and AMDGPU.jl expose more stable
  allocation accounting
- expand backend-native support for map-heavy helper paths where real workloads
  justify it
- keep benchmark reports split by steady-state runtime, cold-start / TTFx,
  precision, and batch throughput
- validate optional self-hosted CUDA/AMDGPU CI runners when those runners become
  part of the project infrastructure

## RTC / HIL DM-To-Pixels Plan
The near-real-time use case is receiving raw deformable-mirror actuator commands
from an RTC and returning detector pixels. The preferred design is to keep
actuator-to-surface reconstruction outside the generic PROPER compatibility
surface when a dedicated AO package can do that job.

Recommended data path:
1. RTC emits raw actuator commands.
2. `AdaptiveOpticsSim.jl` or an application-specific reconstructor converts
   actuator commands into a wavefront-sampled DM surface.
3. `Proper.jl` receives a typed payload through ordinary Julia keywords. Use
   `prop_add_phase(wf, opd_map)` when the AO runtime provides total sampled
   OPD, or `prop_dm(wf, dm_surface)` when the Proper prescription intentionally
   owns the sampled DM surface.
4. Prepared `Proper.jl` execution performs propagation and detector-plane
   sampling with preallocated contexts.

This keeps `Proper.jl` focused on wave-optics propagation and avoids making the
upstream-compatible `prop_dm` actuator/influence-function adapter the real-time
hot path.

Implementation priorities:
1. GPU-native coordinate-grid cubic convolution for arbitrary map projection.
   This is required by map reprojection and the full upstream-compatible
   actuator-space `prop_dm` adapter.
2. GPU-friendly DM surface application remains through `prop_dm(wf, dm_map)`;
   if actuator-space commands must be handled inside `Proper.jl`, route the
   projection through the coordinate-grid interpolation kernel.
3. GPU Zernike map generation only if per-frame Zernike updates are part of a
   real workload.
4. Keep FITS I/O, PSD map generation, fitting routines, plotting, and
   calibration-style helpers host-side unless profiling shows they are inside a
   repeated run loop.

Status:
- Core propagation and prepared batch execution are GPU-capable.
- GPU-native coordinate-grid cubic convolution is implemented for backend
  arrays with explicit failure for unsupported scalar point sampling.
- Full actuator-space `prop_dm` remains a compatibility adapter, not the
  recommended HIL hot path.

The user-facing integration contract is documented in
[`ADAPTIVE_OPTICS_INTEGRATION.md`](ADAPTIVE_OPTICS_INTEGRATION.md).

## Benchmark Commands
```bash
./scripts/benchmark_cpu_gpu.sh
./scripts/profile_core_cpu_gpu.sh
```

For Python-baseline parity plus the full Python-vs-Julia benchmark lane, run:

```bash
PYPROPER_ROOT=/path/to/proper_v3.3.4_python ./scripts/benchmark_all.sh
```
