# Implementation Progress

## Status Snapshot
- Date: 2026-03-13
- Overall: Phase 9 complete; Refactor Track complete; backend extension and CUDA optimization workstream active

## Phase Checklist
- [x] Phase 0: Preflight decisions and contracts accepted
- [x] Phase 1: Foundation and infrastructure
- [x] Phase 2: Core wavefront/math kernels
- [x] Phase 3: I/O and state systems
- [x] Phase 4: Apertures/masks/shapes
- [x] Phase 5: Advanced propagation/DM/maps
- [x] Phase 6: Parity harness and examples
- [x] Phase 7: Benchmarks, CI, release prep
- [x] Phase 8: Full parity closure (Python baseline)
- [x] Phase 9: MATLAB semantic reconciliation and final validation

## Refactor Track Status
- [x] R1: Typed option boundaries (boundary normalization and typed option structs added for propagation/map/geometry families).
- [x] R2: Trait-driven kernel routing (CPU/GPU-ready dispatch wiring, context-routed propagation/interpolation kernels, and optional CUDA smoke test coverage).
- [x] R3: Mutating kernels and workspace reuse across interpolation/geometry hot paths (added `prop_rotate!`, `prop_magnify!`, `prop_cubic_conv_grid!`, geometry `*_!` variants, and reusable interpolation/FFT workspace caches in `RunContext`).
- [x] R4: WaveFront state typing and dispatch simplification (state enums + typed transition selectors for propagation/lens paths; strict typed state transitions).
- [x] R5: Expanded inference/allocation gates and benchmark matrix updates (new R5 gate suite; benchmark matrix expanded for phase-2 kernels, refactor kernel deltas, and end-to-end example workflows; steady-state vs TTFx separation retained).
- [x] GPU extension scaffold: optional `CUDA.jl` package extension with `CuArray` trait registration and public KA-routed interpolation/end-kernel entry points.
- [x] Optional CUDA benchmark lane: separate availability-gated steady-state and supported-kernel reports integrated into `scripts/benchmark_all.sh`.
- [x] Geometry/sampling KA pilot: trait-routed geometry mask kernels plus `prop_szoom!` / `prop_pixellate!`, with CPU pilot benchmarks and CUDA benchmark coverage.
- [x] CUDA hot-path cleanup pass: direct CUDA phase kernels for `prop_qphase` / `prop_ptp`, async KA helper routing, and bounded geometry launches.

## Current Workstream (Phase 1)
- [x] Core policy/trait/context types created.
- [x] One-to-one source file skeleton generated from Python modules.
- [x] One-to-one example file skeleton generated from Python examples.
- [x] Baseline tests added for context and file mapping.
- [ ] Implement `prop_run`/`prop_run_multi` dispatch scaffolding.
- [x] Implement `prop_run`/`prop_run_multi` dispatch scaffolding.
- [x] Add benchmark harness skeleton scripts.

## Current Workstream (Phase 2/3)
- [x] Implemented full Phase-2 propagation kernels and Gaussian-beam state updates:
  - `prop_select_propagator`, `prop_propagate`
  - `prop_ptp`, `prop_stw`, `prop_wts`, `prop_qphase`
  - `prop_lens`
- [x] Added Phase-2 inference checks in tests and dedicated kernel benchmark script.
- [x] Implemented first-pass geometry/mask kernels and obscurations.
- [x] Implemented first-pass map transforms: `prop_magnify`, `prop_rotate`, `prop_resamplemap`.
- [x] Implemented FITS/map baseline wrappers: `prop_fits_read`, `prop_fits_write`, `prop_readmap`, `prop_writemap`, `prop_errormap`.
- [x] Replaced all `_not_implemented` module stubs with callable fallback implementations.

## Current Workstream (Phase 6/7)
- [x] Replaced all example placeholders with runnable Julia scripts (one-to-one filenames).
- [x] Added parity harness skeleton under `test/parity/`.
- [x] Added benchmark scripts for Python steady-state, Julia steady-state, and Julia cold-start.
- [x] Added benchmark summarization script and driver shell script.
- [x] Port example internals to fully match upstream Python physics (all 23 example filenames now contain translated prescription/demo logic).
- [x] Wire benchmark/parity scripts into CI and publish artifacts.

## Next Workstream
- [x] Replace fallback propagation internals with parity implementations.
- [x] Replace fallback PSD/Zernike/polygon/segmented-optics implementations with parity implementations (threshold-gated parity validated against patched Python baseline).
- [x] Replace interpolation and zoom placeholders with cubic-convolution implementations (`libcconv`, `libcconvthread`, `libszoom`, `prop_cubic_conv`).
- [x] Run full example parity report and close threshold gaps against patched Python baseline (combined abs/rel threshold policy and CI gating are green; residual deep-null coronagraph deltas are documented and accepted).
- [x] Document Phase 8 closure evidence (`docs/PHASE8_CLOSURE.md`).
- [x] Perform MATLAB/manual semantic reconciliation on known disagreement hotspots.
- [x] Publish Phase 9 reconciliation report (`docs/PHASE9_RECONCILIATION.md`) and migration guide (`docs/MIGRATION_GUIDE.md`).
- [ ] Complete backend-aware workspace/device cache refactor for CUDA (`docs/CUDA_OPTIMIZATION_PLAN.md`, C2).
- [ ] Reuse CUDA scratch/FFT state in propagation hot paths (`docs/CUDA_OPTIMIZATION_PLAN.md`, C3).
- [ ] Remove remaining host-staged mask/map paths on CUDA (`docs/CUDA_OPTIMIZATION_PLAN.md`, C4).
- [ ] Split CUDA benchmark interpretation by precision regime (`docs/CUDA_OPTIMIZATION_PLAN.md`, C5).

## Latest Pass (2026-03-04)
- [x] Corrected `prop_sinc` to Python-compatible `sin(x)/x`.
- [x] Upgraded `prop_psd_errormap` parity behavior:
  - Python-compatible rotation/inclination quirk handling.
  - Python-compatible `MAX_FREQUENCY` behavior.
  - `FILE` read/reuse and FITS writeout header behavior.
- [x] Added deterministic RNG seeding in parity metric generation scripts.
- [x] Expanded parity reports to include absolute error fields and denominator-floored relative metrics for high-contrast null regions.
- [x] Added executable parity threshold policy (`docs/parity_thresholds.md` + `test/parity/thresholds/example_metrics_thresholds.json`) and CI gating for multi-example parity checks.
- [x] Added hotspot reconciliation tests and completed full-state restore semantics in `prop_state`.
- [x] Added Phase 9 semantic reconciliation report and migration guide.
- [x] Completed R3 refactor scope:
  - introduced `src/core/workspace.jl` and threaded reusable interpolation/FFT scratch through `RunContext`
  - switched hot map kernels to mutating grid/interpolation pathways
  - added mutating geometry mask builders and parity tests for non-mutating wrappers
  - added R3 regression tests for workspace reuse and mutating-kernel parity (`test/test_r3_mutating_workspace.jl`)
- [x] Completed R4 refactor scope:
  - replaced wavefront propagation state `Symbol`s with typed enums in `WaveFront`
  - removed runtime string/symbol transition assembly in `prop_select_propagator`/`prop_lens`
  - routed propagation transition execution through typed selectors in `prop_propagate`
  - added state-typing regression tests (`test/test_r4_state_typing.jl`)
- [x] Completed R5 refactor scope:
  - added expanded inference/allocation gates across propagation, PSD/map, and geometry hotspots (`test/test_r5_performance_gates.jl`)
  - expanded benchmark matrix with dedicated phase-2 kernel report, refactor wrapper-vs-mutating delta report, and end-to-end example workflow report
  - updated benchmark driver and summary reporting while preserving steady-state vs cold-start/TTFx separation
- [x] Added `ProperCUDAExt` package extension:
  - `CUDA.jl` is now a weak dependency rather than a hard package dependency
  - `CuArray` backend/FFT/interpolation traits register through `ext/ProperCUDAExt.jl`
  - public `prop_cubic_conv_grid!`, `prop_rotate!`, `prop_end!`, `prop_qphase`, and mask-application paths avoid host fallbacks when CUDA is available

## Notes
- Runtime compatibility mode flags were removed; parity behavior is anchored to the patched Python baseline (`D-0035`).
- Benchmark policy keeps steady-state comparisons separate from TTFx (`D-0029`).
- CUDA optimization is now tracked separately in `docs/CUDA_OPTIMIZATION_PLAN.md` to keep backend work distinct from the completed parity/refactor milestones.
