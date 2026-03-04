# Implementation Progress

## Status Snapshot
- Date: 2026-03-04
- Overall: In progress

## Phase Checklist
- [x] Phase 0: Preflight decisions and contracts accepted
- [~] Phase 1: Foundation and infrastructure (in progress)
- [ ] Phase 2: Core wavefront/math kernels
- [ ] Phase 3: I/O and state systems
- [ ] Phase 4: Apertures/masks/shapes
- [ ] Phase 5: Advanced propagation/DM/maps
- [ ] Phase 6: Parity harness and examples
- [ ] Phase 7: Benchmarks, CI, release prep

## Current Workstream (Phase 1)
- [x] Core policy/trait/context types created.
- [x] One-to-one source file skeleton generated from Python modules.
- [x] One-to-one example file skeleton generated from Python examples.
- [x] Baseline tests added for context and file mapping.
- [ ] Implement `prop_run`/`prop_run_multi` dispatch scaffolding.
- [ ] Add benchmark harness skeleton scripts.

## Notes
- `compat_mode` is constructor-only and resolves once to a policy type (`D-0017`).
- Benchmark policy keeps steady-state comparisons separate from TTFx (`D-0029`).
