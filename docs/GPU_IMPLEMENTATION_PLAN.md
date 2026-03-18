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

Completed in slice 2:
- the warmed GPU allocation helpers now perform a second warm call before
  measurement so the gates track steady-state behavior rather than first-launch
  backend churn after extension load
- this specifically stabilized the AMDGPU `prop_qphase` gate, whose steady-state
  path is low-kilobyte launch overhead rather than the one-time first-launch
  host churn seen immediately after precompile/extension load

### G4: Promote Mutating, Context-Aware APIs As The GPU Performance Surface
Status: Completed

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

Completed in slice 1:
- public docstrings for `prop_resamplemap!`/`prop_resamplemap` and
  `prop_magnify!`/`prop_magnify` now call out that `ctx` + mutating entry
  points are the repeated GPU performance surface
- user-facing docs now distinguish:
  - warmed steady-state benchmark rows
  - convenience wrappers
  - explicit `RunContext` and prepared execution as the intended GPU contract

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
Status: Completed

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

Completed in slice 1:
- [backend_traits.md](backend_traits.md) now classifies the current backend
  surface into:
  - fully backend-native paths
  - paths that are backend-native only for specific entry points or options
  - intentionally unsupported GPU combinations
- the root README now summarizes the GPU usage contract and links to the full
  support matrix

### G7: Synthetic Core GPU Benchmark And Profiling Harness
Status: Completed

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

Completed in slice 1:
- added a synthetic core propagation-tail workload under
  [`bench/common/core_propagation_tail.jl`](../bench/common/core_propagation_tail.jl)
  built from repeated `prop_lens`/`prop_propagate` transitions
- added steady-state benchmark entry points for:
  - CPU
  - CUDA
  - AMDGPU
- added a profiling driver:
  [`scripts/profile_core_cpu_gpu.sh`](../scripts/profile_core_cpu_gpu.sh)
- integrated the synthetic core workload into the Julia CPU/GPU summary as
  `Synthetic Core Propagation Tail`
- current measured summary row on this machine:
  - Julia CPU: `72.43 ms`
  - Julia AMDGPU: `10.82 ms`
  - CPU/AMDGPU: `6.69x`
- current profile artifacts written to `bench/reports/`:
  - `core_propagation_tail_profile_cpu.txt`
  - `core_propagation_tail_profile_cuda.txt` when CUDA is available
  - `core_propagation_tail_profile_amdgpu.txt`
- current hotspot takeaway:
  - CPU remains dominated by FFT execution in `prop_wts` / `prop_stw`
  - AMDGPU is similarly dominated by rocFFT execution, with `prop_qphase`
    present but secondary

### G8: Julia-Language Cleanup For GPU Paths
Status: In Progress

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

Completed in slice 1:
- replaced broadcast scaling in the warmed FFT propagation path with an
  explicit backend kernel for CUDA/AMDGPU and `rmul!` on the strided CPU path:
  - [`prop_ptp.jl`](../src/prop_ptp.jl)
  - [`prop_wts.jl`](../src/prop_wts.jl)
  - [`prop_stw.jl`](../src/prop_stw.jl)
  - [`ka_kernels.jl`](../src/core/ka_kernels.jl)
- measured effect on this machine after the change:
  - synthetic core propagation tail:
    - AMDGPU host allocations: `1827 -> 1547`
    - AMDGPU host bytes: `62.39 KiB -> 49.58 KiB`
  - supported kernel lane:
  - `prop_ptp`: AMDGPU allocations `514 -> 456`, bytes `14.55 KiB -> 11.80 KiB`
  - `prop_wts`: AMDGPU allocations `433 -> 404`, bytes `11.33 KiB -> 9.95 KiB`
  - `prop_stw`: AMDGPU allocations `433 -> 404`, bytes `11.33 KiB -> 9.95 KiB`

Completed in slice 2:
- directly measured the remaining AMDGPU `prop_qphase` host allocation path and
  found that:
  - first measured call after extension/backend warmup can still be large
  - steady-state repeated `prop_qphase` drops to low-kilobyte host allocation
  - the remaining cost is dominated by kernel launch/runtime overhead rather
    than type instability or a CPU fallback bug
- consequence:
  - `prop_qphase` is no longer the best next GPU optimization target
  - larger returns are more likely in FFT-heavy propagation, geometry/mask
    kernels, or further launch/materialization cleanup around those paths

Completed in slice 3:
- replaced the same-backend full-size GPU `prop_end!` path with the existing
  single-kernel KA copy path instead of the old four-quadrant broadcast/view
  path in [`prop_end.jl`](../src/prop_end.jl)
- measured effect on this machine:
  - direct warmed AMDGPU `prop_end!` allocation:
    - `21.47 KiB -> 8.91 KiB`
  - supported-kernel lane:
    - AMDGPU `prop_end_mutating`:
      - allocs `397 -> 179`
      - bytes `15.58 KiB -> 4.00 KiB`
      - median `132.44 us -> 125.54 us`
- consequence:
  - `prop_end!` is no longer a lagging GPU row; it is now clearly faster than
    the CPU row on this machine

Completed in slice 4:
- replaced the GPU rectangle-mask path in [`ka_kernels.jl`](../src/core/ka_kernels.jl)
  with a single full-grid kernel instead of the old `fill!` plus bounded-box
  kernel sequence
- measured effect on this machine:
  - direct warmed AMDGPU `prop_rectangle!` allocation:
    - `12.91 KiB -> 10.34 KiB`
  - supported-kernel lane:
    - AMDGPU `prop_rectangle_mutating`:
      - allocs `269 -> 199`
      - bytes `7.03 KiB -> 4.86 KiB`
      - median stayed roughly flat (`131.91 us -> 133.12 us`)
- consequence:
  - this was still worth keeping because it reduces host churn materially
    without changing semantics or making the kernel path more complex

Completed in slice 5:
- fused the paired GPU affine-axis fills used by
  [`prop_resamplemap.jl`](../src/prop_resamplemap.jl) and the `QUICK=true`
  path in [`prop_magnify.jl`](../src/prop_magnify.jl) into a single backend
  kernel via new helpers in:
  - [`workspace.jl`](../src/core/workspace.jl)
  - [`ka_kernels.jl`](../src/core/ka_kernels.jl)
- measured effect on this machine:
  - direct warmed AMDGPU `prop_resamplemap!` allocation:
    - `12.14 KiB -> 10.99 KiB`
  - direct warmed AMDGPU quick `prop_magnify!` allocation:
    - `13.89 KiB -> 11.89 KiB`
  - supported-kernel lane:
    - AMDGPU `prop_resamplemap_mutating`:
      - allocs `432 -> 395`
      - bytes `11.08 KiB -> 9.88 KiB`
      - median `439.85 us -> 399.36 us`
    - AMDGPU `prop_magnify_quick_mutating`:
      - allocs `267 -> 199`
      - bytes `7.30 KiB -> 5.56 KiB`
      - median `289.45 us -> 260.97 us`
- consequence:
  - the remaining resample/magnify GPU overhead is no longer dominated by the
    paired affine-axis launches, so the next worthwhile cleanup should move to
    whichever helper still shows the largest steady-state host churn in the
    supported-kernel lane

Completed in slice 6:
- fused the paired GPU damped-sinc table fills used by
  [`prop_szoom.jl`](../src/prop_szoom.jl) into a single backend kernel in
  [`ka_kernels.jl`](../src/core/ka_kernels.jl)
- measured effect on this machine:
  - direct warmed AMDGPU `prop_szoom!` allocation:
    - `12.80 KiB -> 10.80 KiB`
  - direct warmed GPU table-fill split:
    - fused table fill steady-state `8.82 KiB`
    - apply kernel steady-state `8.80 KiB`
  - supported-kernel lane:
    - AMDGPU `prop_szoom_mutating`:
      - allocs `468 -> 427`
      - bytes `12.80 KiB -> 11.23 KiB`
      - median `809.76 us -> 730.13 us`
- consequence:
  - `prop_szoom!` remains one of the heavier non-FFT GPU helpers, but the
    paired table-fill launches are no longer the main avoidable cost

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
- 2026-03-18: `G7` completed with a synthetic core propagation benchmark and
  profiling harness so GPU work no longer depends on WFIRST-derived wrappers.
- 2026-03-18: `G8` slice 1 replaced broadcast scaling in the warmed FFT
  propagation path with explicit backend scaling kernels and reduced AMDGPU
  host allocations in the propagation benchmarks.
- 2026-03-18: `G3` gates were tightened to double-warm before measuring so GPU
  allocation checks reflect the intended steady-state contract.
- 2026-03-18: `G8` slice 5 fused paired GPU affine-axis fills for resample and
  quick magnify, reducing AMDGPU host churn and improving both supported-kernel
  rows.
- 2026-03-18: `G8` slice 6 fused paired GPU damped-sinc table fills for
  `prop_szoom!`, reducing AMDGPU host churn and improving the `prop_szoom`
  supported-kernel row.
- 2026-03-18: direct AMDGPU `prop_qphase` measurement confirmed the remaining
  warm-path host overhead is launch/runtime churn, not a semantics or type
  stability defect.
- 2026-03-18: full-size GPU `prop_end!` switched from quadrant broadcast/view
  staging to a single backend kernel, materially reducing warmed AMDGPU host
  allocations.
- 2026-03-18: GPU rectangle-mask staging collapsed from `fill!` plus bounded
  kernel launch to one full-grid kernel, reducing warmed AMDGPU host
  allocations.
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
- 2026-03-18: G4 completed.
  - public docstrings now make `ctx` + mutating APIs the explicit repeated GPU
    performance surface for `prop_resamplemap` and `prop_magnify`
  - README/docs now state that BenchmarkTools-based GPU benchmark rows are
    warmed steady-state measurements, not cold-start / TTFx
- 2026-03-18: G6 completed.
  - `backend_traits.md` now contains an explicit backend support matrix for
    fully-native, conditionally-native, and intentionally unsupported GPU paths
  - root docs now link that support matrix from the user-facing GPU benchmark
    section
