# Implementation Progress

## Status Snapshot
- Date: 2026-03-04
- Overall: In progress

## Phase Checklist
- [x] Phase 0: Preflight decisions and contracts accepted
- [x] Phase 1: Foundation and infrastructure
- [x] Phase 2: Core wavefront/math kernels
- [~] Phase 3: I/O and state systems (in progress)
- [~] Phase 4: Apertures/masks/shapes (in progress)
- [~] Phase 5: Advanced propagation/DM/maps (in progress)
- [~] Phase 6: Parity harness and examples (in progress)
- [~] Phase 7: Benchmarks, CI, release prep (in progress)
- [ ] Phase 8: Full parity closure (Python baseline)
- [ ] Phase 9: MATLAB semantic reconciliation and final validation

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
- [ ] Port example internals to fully match upstream Python physics.
- [x] Wire benchmark/parity scripts into CI and publish artifacts.

## Next Workstream (Phase 8/9)
- [ ] Replace fallback propagation internals with parity implementations.
- [ ] Replace fallback PSD/Zernike/polygon/segmented-optics implementations with parity implementations.
- [ ] Run full 23-example parity report and close threshold gaps in `:python334`.
- [ ] Perform MATLAB/manual semantic reconciliation on known disagreement hotspots.

## Notes
- `compat_mode` is constructor-only and resolves once to a policy type (`D-0017`).
- Benchmark policy keeps steady-state comparisons separate from TTFx (`D-0029`).
