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
Status: Completed
- Audit remaining keyword-heavy wrappers and move option parsing out of hot kernels.
- Add explicit function barriers where mixed-type inputs currently reach loops.
- Targets:
  - eliminate residual `Any` paths in hot-stack `@code_warntype`.

### O5: Trait-Specialized Fast Paths
Status: Completed
- Convert runtime `isa` branches in hot kernels to method-level dispatch where practical.
- Keep generic fallback methods for non-strided / non-FFTW backends.
- Targets:
  - cleaner IR
  - fewer branch costs in hot loops

### O6: Benchmark and Gate Tightening
Status: Completed
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

### O9: KA Interpolation Pilot
Status: Completed
- Add KA-backed pilot kernels for interpolation-heavy paths:
  - `prop_cubic_conv_grid!`
  - `prop_rotate!` (`linear` and `cubic` methods)
  - inherited public-path coverage for `prop_resamplemap!` and `prop_magnify!(; QUICK=true)` through cubic-grid routing
- Keep loop baselines as the correctness and CPU-performance reference.
- Benchmark KA vs loop directly before enabling any CPU default route.
- Targets:
  - exact parity against current loop implementations
  - explicit benchmark evidence for CPU enable/disable policy
  - trait-routed scaffolding ready for future GPU/backend work

### O10: KA Geometry and Sampling Pilot
Status: Completed
- Add KA-backed pilot kernels for:
  - `prop_rectangle!`
  - `prop_ellipse!`
  - `prop_irregular_polygon!` / `prop_polygon!`
  - `prop_rounded_rectangle!`
  - `prop_szoom!`
  - `prop_pixellate!`
- Keep loop baselines as the correctness and CPU-performance reference.
- Benchmark KA vs loop directly before enabling any CPU default route.
- Targets:
  - exact parity against current loop implementations
  - GPU-ready geometry/sampling kernels without regressing CPU defaults

### O11: CUDA Hot-Path Cleanup Pass
Status: Completed
- Implement the first CUDA-focused runtime cleanup slice:
  - direct KA phase application for `prop_qphase`
  - direct KA frequency-domain phase application for CUDA `prop_ptp`
  - remove unconditional synchronization from low-level KA helpers
  - restrict rectangle / rounded-rectangle KA launches to the affected bounding box
- Keep deeper device-workspace refactors separate and tracked in `docs/CUDA_OPTIMIZATION_PLAN.md`.
- Targets:
  - reduce CUDA launch overhead on small kernels
  - remove CPU-staged map generation from the most obvious propagation hotspots
  - preserve existing CPU behavior and parity results

### O12: Backend-Aware Interp/Mask Workspace Slice
Status: Completed
- Make interpolation-axis workspace buffers backend-aware.
- Make reusable mask buffers backend-aware.
- Update `WaveFront` and `RunContext` constructors to create workspace state from the field/output backend.
- Add typed axis-fill execution routing so CuArray interpolation axes are filled without scalar indexing.
- Keep FFT workspace CPU-owned in this slice; full device-side FFT/cache work stays in `docs/CUDA_OPTIMIZATION_PLAN.md` C3.
- Targets:
  - remove repeated host staging for interpolation axes and aperture mask buffers
  - preserve CPU allocation/performance behavior while enabling backend-preserving scratch on CUDA-visible paths

### O13: Propagation FFT Scratch Reuse Slice
Status: Completed
- Preserve backend complex FFT scratch in `FFTWorkspace`.
- Route CUDA propagation paths for:
  - `prop_ptp`
  - `prop_wts`
  - `prop_stw`
  through workspace-backed in-place `fft!` / `bfft!` transforms.
- Keep CPU FFTW planned path unchanged.
- Targets:
  - reduce temporary transform allocations on CUDA propagation paths
  - keep propagation dispatch type-stable and backend-specific

### O14: Scoped Run-Context Reuse
Status: Completed
- Add task-scoped `RunContext` reuse through `prop_run` so existing prescriptions can reuse backend/workspace state without signature changes.
- Extend `prop_begin` / `prop_wavefront` to honor explicit or scoped workspace/context injection.
- Targets:
  - preserve familiar `prop_run` usage while enabling repeated-run cache reuse
  - keep workspace reuse explicit at the entry-point boundary rather than hidden in hot kernels

### O15: Explicit FFT Planning Policy
Status: Completed
- Make FFT planning an explicit typed policy on `RunContext`.
- Track cached plan flags in `FFTWorkspace` so CPU FFTW plans rebuild when planning mode changes.
- Keep default behavior parity-safe (`FFTEstimateStyle()`), with an opt-in `FFTMeasureStyle()`.
- Targets:
  - make FFT planning behavior visible and testable
  - preserve concrete hot-path dispatch

### O16: Direct Shifted-Ellipse Aperture Path
Status: Completed
- Factor ellipse geometry into a typed descriptor reused by mask-generation and field-application paths.
- Route circular/elliptical aperture and obscuration wrappers through a direct shifted-ellipse path on KA backends.
- Targets:
  - remove redundant mask materialization on CUDA-visible aperture workflows
  - keep wrapper semantics identical to the existing mask path

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
- 2026-03-13: O14 completed.
  - added scoped `RunContext` reuse via `prop_run(...; context=ctx)`.
  - `prop_begin` and `prop_wavefront` now honor explicit/scoped workspace injection and preserve backend/type through reusable workspace state.
- 2026-03-13: O15 completed.
  - `RunContext` now carries typed FFT planning policy (`FFTEstimateStyle` / `FFTMeasureStyle`).
  - `FFTWorkspace` now tracks cached FFTW plan flags and rebuilds plans when policy changes.
- 2026-03-13: O16 completed.
  - circular/elliptical aperture and obscuration wrappers now use a direct shifted-ellipse application path on KA backends.
  - added regression coverage to prove the direct path matches the old mask semantics.
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
- 2026-03-05: O4/O5 completed (first dispatch-focused pass).
  - replaced runtime `isa` branching in `prop_wts`/`prop_stw` with FFT-style + field-type dispatch barriers.
  - replaced runtime strided checks in `prop_end!` with dispatch-routed copy helpers (`noabs` and intensity paths).
  - replaced runtime strided check in PSD map generation with dispatch (`_build_psd_map_default!`).
  - added inference/allocation gates for `prop_wts` and `prop_stw` in `test/test_r5_performance_gates.jl`.
  - post-pass benchmark snapshot (`scripts/benchmark_all.sh`):
    - Python median: `8.4526305e7 ns`
    - Julia median: `2.6214626e7 ns`
    - Python/Julia ratio: `3.224`
- 2026-03-05: O6 completed.
  - benchmark suite rerun and reported through `scripts/benchmark_all.sh` and report summary.
  - performance gates tightened in `test/test_r5_performance_gates.jl` for PSD wrapper/mutating paths and propagation kernels.
- 2026-03-11: O9 completed.
  - added KA-backed internal kernels for cubic-grid interpolation and rotate (`linear`/`cubic`) in `src/core/ka_kernels.jl`.
  - added trait/guard wiring for interpolation kernels in `src/core/traits.jl`; CPU defaults remain disabled (`KA_CUBIC_GRID_MIN_ELEMS = typemax(Int)`, `KA_ROTATE_MIN_ELEMS = typemax(Int)`).
  - added parity tests for KA interpolation kernels vs loop baselines in `test/test_r2_trait_routing.jl`.
  - added dedicated benchmark report `bench/julia/steady_state/ka_interp_kernels.jl` and summary integration.
  - benchmark outcome on CPU (`n = 256`):
    - `cubic_conv_grid`: loop/KA speed ratio `0.870` (KA slower)
    - `rotate_cubic`: loop/KA speed ratio `0.963` (KA slower)
    - `rotate_linear`: loop/KA speed ratio `0.921` (KA slower)
  - decision: keep KA interpolation routes implemented but disabled by default on CPU until CUDA or another backend justifies enabling them.
- 2026-03-12: optional CUDA benchmark lane added.
  - `scripts/benchmark_all.sh` now runs separate CUDA benchmark scripts after the CPU benchmark set.
  - new reports:
    - `bench/reports/julia_cuda_steady_state.json`
    - `bench/reports/cuda_supported_kernels.json`
  - behavior is availability-gated:
    - CUDA available: record synchronized GPU timings for the currently supported GPU subset
    - CUDA unavailable: emit explicit skipped reports instead of failing the benchmark driver
- 2026-03-12: O10 completed for geometry/sampling KA pilot.
  - added trait-routed KA geometry kernels for:
    - `prop_rectangle!`
    - `prop_ellipse!`
    - `prop_irregular_polygon!` / `prop_polygon!`
    - `prop_rounded_rectangle!`
  - added trait-routed sampling kernels for:
    - `prop_szoom!`
    - `prop_pixellate!`
    - `prop_magnify!(; QUICK=false)` now uses `prop_szoom!` instead of allocate-and-copy wrapper flow.
  - all public geometry mask wrappers now preserve the input backend via `similar(wf.field, ...)`.
  - added direct KA-vs-loop equivalence tests in `test/test_r2_trait_routing.jl` and mutating API parity/allocation tests in `test/test_r3_mutating_workspace.jl` / `test/test_r5_performance_gates.jl`.
  - added CPU pilot benchmark report `bench/julia/steady_state/ka_geometry_sampling_kernels.jl`.
  - benchmark outcome on CPU:
    - `rectangle`: loop/KA speed ratio `0.695` (KA slower)
    - `ellipse`: loop/KA speed ratio `0.217` (KA slower)
    - `irregular_polygon`: loop/KA speed ratio `0.988` (near-neutral, KA still slightly slower)
    - `rounded_rectangle`: loop/KA speed ratio `0.516` (KA slower)
    - `szoom`: loop/KA speed ratio `0.507` (KA slower)
    - `pixellate`: loop/KA speed ratio `0.959` (near-neutral, KA still slightly slower)
  - decision: keep geometry/sampling KA routes enabled for CUDA-capable backends and disabled by default on CPU.
- 2026-03-13: O11 completed.
  - added direct CUDA KA kernels for `prop_qphase` and CUDA `prop_ptp` phase application.
  - replaced generic CUDA `prop_qphase` host-staged `rsqr` map generation with direct device-side phase application.
  - replaced CUDA `prop_ptp` frequency-domain broadcast phase application with direct device-side phase application.
  - removed unconditional `AK.synchronize` from low-level KA helpers so synchronization occurs at API/benchmark boundaries.
  - restricted `prop_rectangle!` and `prop_rounded_rectangle!` KA launches to computed bounding boxes rather than the full image.
  - added CUDA smoke coverage for `prop_ptp`.
  - local validation:
    - `julia --project=. test/runtests.jl`: pass
    - `./scripts/benchmark_all.sh`: pass on non-CUDA machine, CUDA lane skipped cleanly
  - follow-up work is tracked in `docs/CUDA_OPTIMIZATION_PLAN.md` for backend-aware workspace/device-cache refactors.
- 2026-03-13: O12 completed.
  - `InterpWorkspace` now preserves the requested backend type instead of always using `Vector`.
  - `MaskWorkspace.mask` now preserves the requested backend type instead of always using `Matrix`.
  - `WaveFront` and `RunContext` now construct these workspaces from the field/output backend.
  - `prop_magnify!` and `prop_resamplemap!` now fill interpolation axes through typed loop/KA dispatch, avoiding host scalar loops for CUDA paths.
  - added CUDA smoke coverage for backend-aware interpolation workspace and mask buffer preservation.
  - local benchmark snapshot:
    - `simple_prescription_256`: `4.01 MiB -> 3.50 MiB`
    - `simple_telescope_256`: `5.02 MiB -> 3.50 MiB`
    - `psdtest_128`: `1.64 MiB -> 1.51 MiB`
- 2026-03-13: O13 completed.
  - `FFTWorkspace` now preserves backend complex scratch buffers while keeping CPU plan caching on the FFTW path.
  - CUDA propagation paths for `prop_ptp`, `prop_wts`, and `prop_stw` now use workspace-backed in-place `fft!` / `bfft!`.
  - added CPU scratch-reuse regression coverage in `test/test_r3_mutating_workspace.jl` and CUDA smoke coverage for backend-preserving FFT scratch in `test/test_r2_trait_routing.jl`.
  - local benchmark snapshot:
    - steady-state CPU: `27.89 ms -> 26.92 ms`
    - phase-2 `prop_ptp`: `12.85 ms -> 10.76 ms`
  - CUDA hardware rerun is still required to measure the actual GPU effect of the in-place transform path.
