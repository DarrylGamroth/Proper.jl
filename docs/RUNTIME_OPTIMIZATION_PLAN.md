# Runtime Optimization Plan

Date: 2026-03-05  
Owner: Proper.jl port effort  
Scope: steady-state runtime/allocations (TTFx explicitly excluded from comparison metrics per D-0029)

## Objectives
- Reduce steady-state allocations and runtime in repeated simulation workloads.
- Preserve executable parity behavior against patched Python 3.3.4 baseline.
- Keep implementation idiomatic Julia: dispatch + typed structs + reusable workspace state.

## Current Baseline (before this plan executes)
- `bench/reports/julia_steady_state.json`: median `3.04447965e7 ns`
- `bench/reports/phase2_kernels.json`:
  - `prop_lens`: `80 B`
  - `prop_qphase`: `0 B`
  - `prop_ptp`: `0 B`
- `bench/reports/refactor_kernels.json`:
  - `psd_errormap_no_apply`: `~0.94 MB`, `~1.01 ms`
- Step-level allocation probe on 512-grid workload:
  - `prop_circular_aperture`: ~4.2 MB
  - `prop_propagate`: ~16.8 MB
  - `prop_end`: ~4.2 MB
  - `prop_wts` / `prop_stw`: ~8.4 MB each

## Design Rules
- Reuse workspace from `WaveFront` for all hot paths when possible.
- Avoid runtime symbol/dict branching in inner loops; normalize options once at boundaries.
- Prefer `*_!` mutating kernels + wrappers.
- Keep parity-first semantics; when behavior changes are needed, update `docs/compat_decisions.md` and tests.

## Phases

### O0: Tracking and Guardrails
Status: Completed
- Add this plan and keep status current.
- Ensure each optimization has:
  - targeted allocation/runtime probe
  - parity check (`test/parity/compare_examples.jl`)
  - regression tests

### O1: FFT-Path Allocation Elimination in Propagation
Status: Completed
- Extend wavefront-owned workspace usage to `prop_wts` and `prop_stw`.
- Use cached FFT plans/scratch buffers (same pattern as `prop_ptp`).
- Targets:
  - `prop_wts`/`prop_stw` hot path: near-zero alloc for strided CPU path.
  - lower `prop_propagate` allocation in representative workloads.

### O2: Mask/Geometry Reuse in Aperture/Obscuration Wrappers
Status: Completed
- Remove repeated full-mask allocations in aperture/obscuration wrappers.
- Reuse a workspace-backed real mask buffer.
- Apply centered mask without `prop_shift_center` temporary arrays.
- Targets:
  - reduce wrapper allocations from MB-scale to KB/near-zero in warm path.

### O3: End-of-Run Output Buffer Reuse
Status: Completed
- Add mutating extraction/copy path for `prop_end` internals (`prop_end!`-style buffer option).
- Keep existing API unchanged; wrapper allocates only when caller does not provide output buffer.
- Targets:
  - reduce repeated `prop_end` allocations in batched workflows.

### O4: Boundary Normalization and Function-Barrier Sweep
Status: Pending
- Audit remaining keyword-heavy wrappers and move option parsing out of hot kernels.
- Add explicit function barriers where mixed-type inputs currently reach loops.
- Targets:
  - eliminate residual `Any` paths in hot-stack `@code_warntype`.

### O5: Trait-Specialized Fast Paths
Status: Pending
- Convert runtime `isa` branches in hot kernels to method-level dispatch where practical.
- Keep generic fallback methods for non-strided / non-FFTW backends.
- Targets:
  - cleaner IR
  - fewer branch costs in hot loops

### O6: Benchmark and Gate Tightening
Status: Pending
- Re-run `scripts/benchmark_all.sh`.
- Update performance gates in tests to reflect improved steady-state allocations.
- Record results and remaining hotspots in this file.

### O7: KA/AK Pilot for Shift Kernels
Status: Completed
- Introduce trait-routed pilot kernels using `KernelAbstractions.jl` with `AcceleratedKernels.jl` backend/sync helpers.
- Scope:
  - shifted mask apply path (`_apply_shifted_mask!`)
  - shifted output copy/intensity path in `prop_end!`
- Keep loop fallbacks as default for smaller arrays via threshold guard.
- Targets:
  - maintain parity and tests
  - validate neutral-or-better runtime on representative workloads
  - establish backend-dispatch scaffolding for future GPU methods

### O8: PSD Map Mutating Path and Benchmark De-biasing
Status: Completed
- Add preallocated-output `prop_psd_errormap!` API and internal output-buffer overloads.
- Keep parity behavior unchanged (same returned map values and application semantics).
- Replace `prop_shift_center` allocation in PSD output normalization with workspace-backed in-place shift.
- Update refactor kernel benchmark to:
  - include PSD wrapper vs mutating pair
  - avoid per-call RNG-object construction noise in PSD timings
- Targets:
  - mutating PSD path near-zero extra allocations (beyond caller-owned output)
  - preserve parity and regression tests

## Execution Log
- 2026-03-05: Plan created. Starting O1.
- 2026-03-05: O1 completed.
  - `prop_wts` warm-path allocation: `0 B` (strided CPU path).
  - `prop_stw` warm-path allocation: `0 B` when staying on spherical path.
  - representative one-step `prop_propagate` allocation probe reduced from ~16.8 MB to ~8.4 MB.
- 2026-03-05: O2 started.
  - added wavefront-owned mask buffer reuse.
  - `prop_circular_aperture` warm-path allocation reduced to ~`1328 B` on 512-grid probe.
- 2026-03-05: benchmark update after O1 + O2 start (`scripts/benchmark_all.sh`).
  - steady-state median improved to `2.6180651e7 ns` (from `3.04447965e7 ns` pre-plan snapshot).
  - Python/Julia ratio moved to `0.488` in this benchmark setup.
  - example workflow allocations reduced significantly:
    - `simple_prescription_256`: `~7.35 MB -> ~4.73 MB`
    - `psdtest_128`: `~4.35 MB -> ~2.38 MB`
    - `simple_telescope_256`: `~12.61 MB -> ~5.79 MB`
- 2026-03-05: O2 completed.
  - wrappers updated to reuse wavefront mask workspace:
    - `prop_circular_aperture` warm-path allocation: `1328 B` (512-grid probe)
    - `prop_elliptical_aperture` warm-path allocation: `1328 B` (512-grid probe)
    - `prop_rectangular_aperture` warm-path allocation: `832 B` (512-grid probe)
- 2026-03-05: O3 completed.
  - added `prop_end!` mutating output path and `prop_end(wf, out; ...)` buffer overload.
  - `prop_end!` warm-path allocation: `0 B` (full and extract paths on 512-grid probe).
  - default allocating wrapper now allocates output-only footprint (`~2.10 MB` for 512x512 intensity output).
- 2026-03-05: benchmark update after O2 + O3 completion (`scripts/benchmark_all.sh`).
  - steady-state median: `2.79100035e7 ns`; Python/Julia ratio in this setup: `0.467`.
  - note: this ratio was from the old non-parity Python toy harness and is superseded by the corrected parity benchmark entry below.
  - example workflow medians / allocations:
    - `simple_prescription_256`: `7.00e6 ns`, `4,206,240 B`
    - `psdtest_128`: `4.05e6 ns`, `2,248,328 B`
    - `simple_telescope_256`: `1.60e7 ns`, `5,264,464 B`
- 2026-03-05: benchmark harness corrected to use Python PROPER workload parity (same 512-grid propagation sequence as Julia benchmark).
  - `bench/python/run.py` now imports `../proper_v3.3.4_python` and executes `prop_begin -> prop_circular_aperture -> prop_lens -> prop_propagate -> prop_end`.
  - `scripts/benchmark_all.sh` now prefers `.venv-parity/bin/python` for baseline execution.
  - updated steady-state comparison after harness fix:
    - Python median: `8.6305942e7 ns`
    - Julia median: `2.7026525e7 ns`
    - Python/Julia ratio: `3.193`
- 2026-03-05: allocation reduction pass for PSD + geometry kernels completed.
  - `prop_psd_errormap` now reuses FFT workspace scratch/plans and phase matrix, and avoids shift-copy temporaries in the generated-map path.
  - `prop_irregular_polygon`/`prop_polygon` now avoid `coordinate_axis`/`collect`/broadcast temporaries in hot loops and reuse workspace vertex buffers.
  - `prop_ellipse` now removes temporary corner/offset vectors in favor of scalar loop math.
  - updated kernel allocation metrics (`bench/julia/steady_state/refactor_kernels.jl`):
    - `psd_errormap_no_apply`: `~941,344 B -> ~154,224 B`
    - `ellipse` mutating: `~1,328 B -> ~832 B`
    - `polygon` mutating: `~6,400 B -> ~832 B`
    - `irregular_polygon` mutating: `~5,408 B -> ~544 B`
- 2026-03-05: O7 KA/AK pilot implemented for shifted kernels.
  - added `ShiftKernelStyle` trait with per-kernel guarded enable (`KA_MASK_MIN_ELEMS`, `KA_END_MIN_ELEMS`).
  - added KA kernel implementations for shifted mask apply and shifted `prop_end!` copy/intensity.
  - integrated trait routing in `prop_circular_aperture` mask application and `prop_end!` strided fast paths.
  - updated CPU defaults after benchmarking: mask path remains loop-default; `prop_end!` shifted copy/intensity keeps KA pilot for large arrays.
- 2026-03-05: O8 completed.
  - added `prop_psd_errormap!(out, wf, ...)` mutating API and internal output-buffer overload.
  - replaced PSD `prop_shift_center` allocation with workspace-backed in-place shift for matrix fast path.
  - tightened PSD performance gates:
    - wrapper allocation threshold: `< 100_000` bytes (64-grid test)
    - mutating allocation threshold: `< 10_000` bytes (64-grid test)
  - updated refactor benchmark report schema to include PSD wrapper/mutating pair + hotspot entries.
