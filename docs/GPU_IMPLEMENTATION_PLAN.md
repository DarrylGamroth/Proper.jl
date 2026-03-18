# GPU Implementation Plan

Working implementation plan for GPU backend quality, performance, and API
honesty across CUDA and AMDGPU. This is the active maintainer plan for
GPU-focused work; backend-specific notes such as
[`CUDA_OPTIMIZATION_PLAN.md`](CUDA_OPTIMIZATION_PLAN.md) remain useful as
historical or backend-local follow-up records.

Date: 2026-03-17  
Owner: Proper.jl port effort  
Scope: shared GPU execution model, propagation-path allocation discipline,
backend-native workspaces, and removal of hidden host fallback in GPU-visible
paths

## Problem Statement
- CUDA and AMDGPU are both functional and useful for propagation-heavy
  workloads.
- The core propagation dispatch is in reasonable shape, but the overall GPU
  surface is still uneven.
- Several GPU-visible paths silently materialize on the CPU or allocate fresh
  context/workspace state in convenience wrappers.
- The next work should improve the shared GPU implementation, not tune
  WFIRST-specific code.

## Goals
- Preserve parity behavior against the patched Python 3.3.4 baseline.
- Keep the core backend-generic through traits and typed contexts.
- Make warmed propagation paths allocation-free or near-allocation-free.
- Remove silent host fallback from GPU-visible hot paths.
- Make the performance contract explicit:
  - mutating `!` APIs + explicit `RunContext` are the performance surface
  - allocating wrappers remain convenience APIs

## Non-Goals
- Do not force zero allocations in one-shot convenience wrappers.
- Do not optimize reference models like WFIRST as model-specific special cases.
- Do not add runtime compatibility flags.
- Do not introduce speculative shared-memory/tiling kernels without profile
  evidence.

## Current Findings Driving This Plan

### F-001: Hidden host fallback remains in GPU-visible APIs
Severity: High

Examples:
- [prop_rotate.jl](../src/prop_rotate.jl)
- [prop_cubic_conv.jl](../src/prop_cubic_conv.jl)
- [prop_end.jl](../src/prop_end.jl)

Current behavior:
- unsupported or generic paths allocate host buffers and call `Matrix(...)`
- results are copied back to the device/backend array

Why it matters:
- this hides performance cliffs
- it makes GPU benchmark results harder to interpret
- it violates the intended backend-preserving execution model

### F-002: `FFTWorkspace` is not fully backend-native
Severity: High

Examples:
- [workspace.jl](../src/core/workspace.jl)

Current behavior:
- `rho2` is stored as `Matrix{T}` even when the rest of the workspace is device
  resident
- generic propagation must `backend_adapt(...)` the map each call

Why it matters:
- it keeps part of the propagation model CPU-owned
- it complicates truly allocation-stable generic backend fallback

### F-003: Convenience wrappers recreate context/workspace state too easily
Severity: Medium

Examples:
- [prop_resamplemap.jl](../src/prop_resamplemap.jl)
- [prop_magnify.jl](../src/prop_magnify.jl)
- implicit `RunContext(wf)` wrappers in propagation/lens helpers

Why it matters:
- GPU users can accidentally benchmark setup churn instead of steady-state work
- allocation behavior looks worse than the real prepared path

### F-004: Map-application paths are still CPU-leaning
Severity: Medium

Examples:
- [prop_add_phase.jl](../src/prop_add_phase.jl)
- [prop_multiply.jl](../src/prop_multiply.jl)
- [prop_divide.jl](../src/prop_divide.jl)
- [prop_errormap.jl](../src/prop_errormap.jl)
- [prop_readmap.jl](../src/prop_readmap.jl)
- [prop_psd_errormap.jl](../src/prop_psd_errormap.jl)

Why it matters:
- many of these paths still allocate shifted temporaries or adapt buffers
- GPU support is functional, but not yet as disciplined as the propagation core

### F-005: GPU regression coverage is benchmark-first, allocation-second
Severity: Medium

Examples:
- [test_r5_performance_gates.jl](../test/test_r5_performance_gates.jl)
- [test_r2_trait_routing.jl](../test/test_r2_trait_routing.jl)
- AMD/CUDA benchmark lanes under `bench/julia/`

Why it matters:
- time-only GPU checks miss warm-path regressions in host allocations,
  context churn, or hidden fallback

## Performance Contract

### What should be zero-allocation
- Warmed `prop_qphase`
- Warmed `prop_ptp`
- Warmed `prop_wts`
- Warmed `prop_stw`
- Warmed `prop_end!`
- Warmed mutating GPU map-application kernels once they are promoted to the
  performance surface

### What does not need to be zero-allocation
- allocating convenience wrappers
- first-call compilation
- first plan creation
- FITS I/O
- model setup/preparation

### Operational definition
- “Zero allocation” means zero steady-state host allocations in the hot loop,
  plus no repeated device-side buffer/plan creation after warmup.

## Workstreams

### G0: Tracking, Scope, and Benchmark Discipline
Status: In Progress

Tasks:
- [ ] Treat this file as the active GPU tracker.
- [ ] Keep backend-specific notes as secondary documents.
- [ ] Separate:
  - correctness smoke
  - steady-state benchmark
  - hotspot profiling
  - cold-start/TTFx
- [ ] Add explicit language to benchmark docs that `BenchmarkTools` covers
  steady-state, not cold-start.

Acceptance:
- Docs consistently distinguish steady-state from cold-start.
- GPU benchmark outputs remain comparable across CPU/CUDA/AMDGPU.

### G1: Remove Hidden Host Fallback From GPU-Visible Hot Paths
Status: Completed

Priority: Highest

Targets:
- [ ] [prop_rotate.jl](../src/prop_rotate.jl)
- [ ] [prop_cubic_conv.jl](../src/prop_cubic_conv.jl)
- [ ] generic [prop_end.jl](../src/prop_end.jl) paths

Rules:
- Prefer true backend implementations.
- If an operation is not supported on GPU, fail explicitly with a clear error.
- Do not silently bounce through `Matrix(...)` in a hot path.

Acceptance:
- No hidden host materialization in warmed GPU benchmark kernels.
- Trait routing clearly distinguishes:
  - native backend execution
  - unsupported backend execution

Completed in slice 1:
- `prop_rotate` now throws on unsupported non-CPU backend/layout combinations
  instead of silently materializing on the host
- `prop_cubic_conv` now throws on unsupported non-CPU point/grid paths instead
  of silently materializing on the host
- `prop_end!` now requires matching output and wavefront backends so host/device
  transfer stays explicit at the call site

### G2: Make Workspaces Fully Backend-Native
Status: Completed

Priority: Highest

Targets:
- [ ] Parameterize `FFTWorkspace.rho2` by backend array type
- [ ] Fill cached frequency maps on the active backend where applicable
- [ ] Audit remaining CPU-owned fields in backend-visible workspace state

Acceptance:
- No CPU-only cache fields remain in GPU propagation workspaces unless the path
  is explicitly CPU-only.
- Generic propagation fallback does not repeatedly adapt CPU maps to the device.

Completed in slice 1:
- `FFTWorkspace.rho2` is now backend-parametric instead of hard-coded as
  `Matrix{T}`
- CUDA and AMDGPU now construct backend-native `rho2` caches
- `ensure_rho2_map!` fills the cached map on the active backend
- optional backend smoke now verifies that `ensure_rho2_map!` returns
  `CuArray`/`ROCArray` storage on those backends

### G3: Establish Warmed Propagation Allocation Gates For GPU
Status: Completed

Priority: High

Targets:
- [ ] Add CUDA/AMDGPU warm-path checks for:
  - `prop_qphase`
  - `prop_ptp`
  - `prop_wts`
  - `prop_stw`
  - `prop_end!`
- [ ] Track host allocations and state churn across repeated warmed calls
- [ ] Add device benchmark notes where direct device allocation accounting is
  not yet practical

Acceptance:
- GPU propagation regression tests fail on repeated warm-path host allocation
  regressions.
- Benchmarks and tests use the same warmed execution assumptions.

Completed in slice 1:
- optional CUDA/AMDGPU smoke now records warmed host-allocation bounds for:
  - `prop_qphase`
  - `prop_ptp`
  - `prop_wts`
  - `prop_stw`
  - `prop_end!`
- the same smoke verifies state-churn invariants that matter for warm execution:
  - FFT plans are reused across repeated propagation calls
  - cached `rho2` storage is reused by `ensure_rho2_map!`
- current gates are bounded host-allocation checks, not a false zero-allocation
  claim for every GPU backend; device allocation tracking remains a benchmark
  and profiling concern for now

### G4: Promote Mutating, Context-Aware APIs As The GPU Performance Surface
Status: Not Started

Priority: High

Targets:
- [ ] Audit convenience wrappers that create `RunContext(typeof(out))` or fresh
  scratch state
- [ ] Ensure docs call out that `ctx` + `!` methods are the performance path
- [ ] Consider lightweight warnings or documentation-only guidance where wrapper
  misuse is common

Examples:
- [prop_resamplemap.jl](../src/prop_resamplemap.jl)
- [prop_magnify.jl](../src/prop_magnify.jl)

Acceptance:
- User-facing docs clearly separate convenience APIs from performance APIs.
- Benchmarks use explicit contexts/prepared execution consistently.

### G5: Convert Map-Application Helpers To Scratch-Backed GPU Paths
Status: Completed

Priority: High

Targets:
- [ ] [prop_add_phase.jl](../src/prop_add_phase.jl)
- [ ] [prop_multiply.jl](../src/prop_multiply.jl)
- [ ] [prop_divide.jl](../src/prop_divide.jl)
- [ ] [prop_errormap.jl](../src/prop_errormap.jl)
- [ ] [prop_readmap.jl](../src/prop_readmap.jl)
- [ ] [prop_psd_errormap.jl](../src/prop_psd_errormap.jl)

Approach:
- prefer mutating centered-shift/application helpers
- reuse backend scratch owned by `ProperWorkspace`
- avoid repeated `backend_adapt(...)` and `copy(...)` in warmed loops

Acceptance:
- GPU-visible map-application paths stay on backend for their hot transforms.
- No silent CPU shift/magnify/rotate staging remains in promoted fast paths.

Completed in slice 1:
- `prop_shift_center!` now has a backend-native KA shift path for CUDA and
  AMDGPU arrays, so same-backend map centering no longer falls back to scalar
  host iteration
- `prop_add_phase`, `prop_multiply`, and `prop_divide` now reuse backend-native
  FFT scratch via `shift_center_for_wavefront!` instead of allocating shifted
  temporaries for same-backend GPU maps
- `prop_readmap` now resamples onto the wavefront backend and returns a
  backend-native shifted map for GPU wavefronts
- `prop_errormap` now keeps rotate/magnify and final application on the
  wavefront backend when the map has already been promoted there
- `prop_psd_errormap` now returns backend-native outputs for non-CPU wavefronts
  and applies shifted maps on-backend instead of bouncing through host-adapted
  temporaries

### G6: Honest Backend Support Matrix
Status: Not Started

Priority: Medium

Tasks:
- [ ] Audit which public operations are:
  - fully backend-native
  - backend-native only for some methods/options
  - currently unsupported on GPU
- [ ] Document that matrix in:
  - [backend_traits.md](backend_traits.md)
  - README-facing benchmark/backend sections
- [ ] Make unsupported combinations fail clearly rather than fallback silently

Acceptance:
- A contributor or user can tell which GPU paths are real.
- Trait tests cover the intended support matrix.

### G7: Synthetic Core GPU Benchmark And Profiling Harness
Status: Not Started

Priority: Medium

Tasks:
- [ ] Add a small core-only propagation benchmark that mirrors the dominant
  propagation sequence identified by WFIRST, without using WFIRST itself as the
  optimization target
- [ ] Keep WFIRST as a validation workload, not the optimization harness
- [ ] Record hotspot snapshots for:
  - CPU
  - CUDA
  - AMDGPU

Acceptance:
- GPU performance work can be evaluated on a shared core workload rather than
  model-specific wrappers.

### G8: Julia-Language Cleanup For GPU Paths
Status: Not Started

Priority: Medium

Tasks:
- [ ] Continue replacing avoidable eager materialization in hot helper paths
- [ ] Keep hot structs fully concrete/parametric
- [ ] Move backend selection entirely into dispatch where runtime branches still
  remain in hot code
- [ ] Prefer reusable buffers over `backend_adapt(...)` in warmed loops

Acceptance:
- No new abstract-field or ad hoc dynamic-dispatch regressions in GPU hot
  paths.
- Hot GPU helpers remain inference-friendly.

## Immediate Execution Order

This is the intended implementation order for the next GPU-focused slices.

1. G1: remove hidden host fallback in `prop_rotate`, `prop_cubic_conv`, and
   generic `prop_end!`
2. G2: make `FFTWorkspace` fully backend-native
3. G3: add warmed GPU propagation allocation/regression gates
4. G5: convert map-application helpers to scratch-backed GPU paths
5. G6-G8: docs/support-matrix/synthetic benchmark cleanup

## Verification Matrix

Every GPU-focused change should include the applicable checks below.

### Required
- [ ] `julia --project=. test/runtests.jl`
- [ ] optional GPU smoke in [test_r2_trait_routing.jl](../test/test_r2_trait_routing.jl)
- [ ] CPU/GPU benchmark lane:
  - `./scripts/benchmark_cpu_gpu.sh`

### When touching propagation
- [ ] warmed repeated benchmark for `prop_qphase`, `prop_ptp`, `prop_wts`,
  `prop_stw`
- [ ] no new host fallback in the changed path

### When touching map/interp paths
- [ ] mutating benchmark lane if applicable
- [ ] parity smoke against CPU reference on representative sizes

### When changing semantics or support matrix
- [ ] update [compat_decisions.md](compat_decisions.md) if behavior changes
- [ ] update [backend_traits.md](backend_traits.md)

## Open Questions
- Should unsupported GPU combinations throw immediately or be documented as
  CPU-only convenience paths? Default answer: throw in hot paths, keep fallback
  only for explicitly convenience-oriented wrappers.
- How much device-allocation accounting do we want in tests? Default answer:
  host-allocation regression checks plus warmed benchmark policy are sufficient
  for now.
- Should we introduce an internal `AbstractFFTs`-style backend layer later?
  Default answer: yes as future work when GPU FFT reuse becomes a portability
  blocker, not before.

## Execution Log
- 2026-03-17: Plan created from GPU implementation review.
- 2026-03-17: Initial findings captured:
  - hidden host fallback in rotate/cubic-conv/end
  - CPU-owned `rho2` cache in `FFTWorkspace`
  - convenience-wrapper context churn in resample/magnify paths
  - GPU map-application paths still rely on temporary shifts/adaptation
- 2026-03-17: G1 slice 1 implemented.
  - removed hidden host fallback from unsupported non-CPU `prop_rotate` paths
  - removed hidden host fallback from unsupported non-CPU `prop_cubic_conv`
    scalar/pointwise/coordinate-grid paths
  - made `prop_end!` require matching output and wavefront backends
  - added regression coverage in `test/test_r2_trait_routing.jl`
- 2026-03-17: G2 completed.
  - `FFTWorkspace.rho2` is now backend-native and no longer hard-coded as a CPU
    `Matrix`
  - CUDA and AMDGPU extensions now construct backend-native cached frequency
  maps alongside backend-native FFT scratch
  - `ensure_rho2_map!` now fills the cache on the active backend instead of
    forcing a CPU fill path
- 2026-03-17: G3 completed.
  - optional backend smoke now includes warmed propagation host-allocation
    regression gates for `prop_qphase`, `prop_ptp`, `prop_wts`, `prop_stw`,
    and `prop_end!`
  - the same smoke now checks repeated warm-path reuse of cached FFT plans and
    `rho2` storage
  - current GPU warm-path checks are bounded-allocation guards rather than
    zero-allocation guarantees, reflecting the measured AMDGPU host-side
    overhead today
- 2026-03-17: G5 completed.
  - added backend-native KA center-shift support for CUDA and AMDGPU arrays
  - `prop_add_phase`, `prop_multiply`, and `prop_divide` now reuse
    backend-native scratch for same-backend map application
  - `prop_readmap`, `prop_errormap`, and `prop_psd_errormap` now keep promoted
    GPU maps on-backend through their hot transform/apply path
  - optional GPU smoke now exercises map application, FITS map reads, error-map
    application, and PSD error-map output/application on AMDGPU
