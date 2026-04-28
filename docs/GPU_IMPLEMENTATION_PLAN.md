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

## Benchmark Commands
```bash
./scripts/benchmark_cpu_gpu.sh
./scripts/profile_core_cpu_gpu.sh
```

For Python-baseline parity plus the full Python-vs-Julia benchmark lane, run:

```bash
PYPROPER_ROOT=/path/to/proper_v3.3.4_python ./scripts/benchmark_all.sh
```
