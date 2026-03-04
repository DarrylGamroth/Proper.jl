# Implementation Progress

## Status Snapshot
- Date: 2026-03-04
- Overall: In progress

## Phase Checklist
- [x] Phase 0: Preflight decisions and contracts accepted
- [~] Phase 1: Foundation and infrastructure (in progress)
- [~] Phase 2: Core wavefront/math kernels (in progress)
- [ ] Phase 3: I/O and state systems
- [ ] Phase 4: Apertures/masks/shapes
- [ ] Phase 5: Advanced propagation/DM/maps
- [~] Phase 6: Parity harness and examples (in progress)
- [~] Phase 7: Benchmarks, CI, release prep (in progress)

## Current Workstream (Phase 1)
- [x] Core policy/trait/context types created.
- [x] One-to-one source file skeleton generated from Python modules.
- [x] One-to-one example file skeleton generated from Python examples.
- [x] Baseline tests added for context and file mapping.
- [ ] Implement `prop_run`/`prop_run_multi` dispatch scaffolding.
- [x] Implement `prop_run`/`prop_run_multi` dispatch scaffolding.
- [x] Add benchmark harness skeleton scripts.

## Current Workstream (Phase 2/3)
- [x] Implemented first-pass optics kernels: `prop_qphase`, `prop_lens`, `prop_propagate`, `prop_ptp`, `prop_wts`, `prop_stw`.
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
- [ ] Wire benchmark/parity scripts into CI and publish artifacts.

## Notes
- `compat_mode` is constructor-only and resolves once to a policy type (`D-0017`).
- Benchmark policy keeps steady-state comparisons separate from TTFx (`D-0029`).
