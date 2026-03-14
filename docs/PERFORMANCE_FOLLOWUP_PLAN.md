# Performance Follow-Up Plan

Date: 2026-03-13  
Owner: Proper.jl port effort  
Scope: resolve the remaining GPU/CPU design questions, implement the immediate runtime fixes, and define the next optimization gates

## Problem Statement
- CUDA is now materially useful for representative workloads, but a few kernel families still underperform CPU or show precision-sensitive behavior.
- The remaining questions are no longer about parity correctness; they are about execution shape, cache lifetime, benchmark interpretation, and how far the compatibility surface should constrain the core.
- The immediate work should fix the highest-value issues without overfitting the implementation to one GPU.

## Guiding Decisions
- Keep an idiomatic Julia core and preserve the familiar `prop_*` API as a compatibility wrapper layer.
- Prefer dispatch and typed runtime context (`RunContext`, backend traits, FFT planning traits) over runtime branching.
- Do not add tiling/shared-memory kernels speculatively. Require profile evidence first.
- Specialize only for stable algorithm families or execution strategies, not for sizes or user-level knobs.
- Keep parity-first behavior against the patched Python 3.3.4 baseline.

## Findings Driving This Plan
1. Type instability is not the main limiter in current hot paths.
2. The largest remaining GPU costs are algorithm-shape and precision-regime issues, not host/device copies in the core propagation path.
3. The public compatibility API is acceptable for users, but the core should continue moving toward typed contexts, reusable workspaces, and prepared execution state.
4. The current circular-aperture CUDA path still pays generic ellipse costs.
5. The existing CUDA precision split mixed FP64 and FP32 steady-state workloads inside one process, which makes interpretation weaker than it should be.

## Immediate Work Items

### F1: Thread Context Through Remaining Hot Internal Callers
Status: Completed
- Add `RunContext`-aware `prop_lens` overloads.
- Keep the existing public API unchanged.
- Make `prop_lens` call `prop_qphase(..., ctx)` instead of constructing a fresh default path.
- Update CUDA steady-state workloads to use the context-aware lens path.

Acceptance:
- Existing `prop_lens(wf, lens_fl[, surface_name])` calls still work.
- `prop_lens(wf, lens_fl, ctx[, surface_name])` and `prop_lens(wf, lens_fl, surface_name, ctx)` are supported.
- Parity and regression tests stay green.

### F2: Reconcile CUDA Steady-State Benchmark Methodology
Status: Completed
- Split CUDA steady-state FP64 and FP32 workload timings into standalone scripts.
- Keep the existing CUDA steady-state report as the default compatibility-facing row.
- Remove mixed steady-state workload timing from the general precision-split script.
- Update skip handling, summary tables, Markdown, and CSV reports.

Acceptance:
- FP64 and FP32 workload rows come from separate Julia invocations.
- `CUDA Precision Split` uses standalone workload reports rather than mixed-process measurements.
- Standard CUDA steady-state output remains available for CPU/Python/CUDA comparison.

### F3: Specialize Circular Aperture/Obscuration GPU Execution
Status: Completed
- Add a direct circle-specific shifted-field kernel for KA backends.
- Keep the generic ellipse implementation as the fallback path.
- Route `prop_circular_aperture` and `prop_circular_obscuration` through the circle path on KA-capable backends.
- Preserve parity semantics, including `NORM`, `DARK`, and obscuration inversion behavior.

Design notes:
- Use FFT-order coordinates directly on the field instead of materializing or sampling a generic ellipse mask.
- Keep the algorithm branch-free where possible for clearly inside/outside pixels.
- Do not yet introduce shared memory or multi-kernel tiling; first get the circle-specific direct path in place and measure it.

Acceptance:
- Exact or threshold-tight parity versus the existing generic circular behavior.
- CUDA smoke coverage exercises the new path with `CUDA.allowscalar(false)`.
- CPU behavior remains unchanged unless the new direct path is explicitly chosen by dispatch.

### F4: Add Profiling Gate For Interpolation Tiling / Shared Memory Work
Status: Completed
- Add a dedicated CUDA interpolation profiling harness.
- Cover:
  - `prop_cubic_conv_grid!`
  - `prop_rotate!`
  - `prop_resamplemap!`
  - `prop_magnify!(; QUICK=true)`
- Capture host-side `Profile.@profile` evidence and device-side `CUDA.@profile` markers on representative workloads.
- Document that tiling/shared-memory follow-up work is gated on those results.

Acceptance:
- A reproducible profiling script exists in-tree.
- The plan explicitly records which kernels are candidates for tiling and which are not.

## Design Review Questions And Default Answers

### Q1: Should we use tiling / `@localmem`?
Default answer: Not yet.
- Use tiling only for kernels with real source-tile reuse.
- Current best candidates are interpolation-family kernels, not FFT-driven propagation.
- Geometry masks and `prop_ptp`/`prop_wts`/`prop_stw` should not be force-fit into shared-memory patterns.

### Q2: Are we iterating in the right order?
Default answer: Mostly yes on CPU; GPU needs profiling rather than intuition.
- CPU loops should keep the first index inner for column-major arrays.
- GPU work should be judged by coalescing and total kernel shape, not by CPU row/column heuristics alone.
- Add a targeted audit item if profiling shows memory throughput bottlenecks.

### Q3: Should kernels be further decomposed?
Default answer: Selectively.
- Prefer specialized common-case kernels over decomposing everything into smaller generic pieces.
- Circular aperture is the immediate example.
- Keep generic ellipse/polygon/etc. as semantic fallbacks.

### Q4: Are we using types, selectors, and traits correctly?
Default answer: Largely yes, but consistency still matters.
- Thread `RunContext` through hot internal callers.
- Keep backend/FFT/interpolation/planning in traits and typed context.
- Do not push dynamic sizes or mutable user options into type parameters.

### Q5: Are we overspecializing?
Default answer: Not materially, but keep guarding against it.
- Specialize on algorithm families and execution backends.
- Avoid value-based type explosion for sizes, flags, or frequently varying options.

### Q6: Is the PROPER API holding the core back?
Default answer: Somewhat, so the core should keep getting cleaner.
- Keep the user-facing `prop_*` compatibility layer.
- Continue moving the implementation core toward reusable contexts, typed options, mutating kernels, and prepared execution state.
- Long-term performance work should target the core first and let the wrappers stay thin.

## Follow-On Items After F1-F4

### N1: Complete context threading audit for remaining hot internal helpers
- Confirm no remaining internal hot-path calls recreate default context unnecessarily.

### N2: Evaluate a specialized `prop_end!` fast path only if it appears in real workflow profiles
Status: In Progress
- Common full-frame `prop_end!` now uses direct quadrant copy/broadcast paths instead of the generic shifted-copy kernel.
- Remaining question is CUDA-specific performance on hardware after this fast path lands.

### N3: Evaluate interpolation tiling only if F4 profiling shows global-memory reuse pressure
Status: Gated
- Likely candidates: cubic convolution / rotate / resample.
- Non-candidates: FFT propagation.

### N4: Define a prepared-simulation core API on top of `RunContext` and `ProperWorkspace`
- Goal: idiomatic Julia core with thin compatibility wrappers.

## Tracking
- This plan supplements `docs/CUDA_OPTIMIZATION_PLAN.md` and `docs/RUNTIME_OPTIMIZATION_PLAN.md`.
- Immediate work in this pass: F1-F4.

## Execution Log
- 2026-03-13: Plan created.
- 2026-03-13: F1 completed.
  - `prop_lens` now has context-aware overloads and threads `ctx` into `prop_qphase`.
  - regression coverage added for context-routed lens behavior.
- 2026-03-13: F2 completed.
  - CUDA steady-state FP64/FP32 workload timing now comes from standalone scripts.
  - the mixed precision-split script now benchmarks kernels only.
  - skip handling and summary reporting now include the standalone workload reports.
- 2026-03-13: F3 completed.
  - circular aperture/obscuration now use a circle-specific KA execution path on KA backends.
  - the generic ellipse implementation remains the fallback path.
  - parity regression coverage added for the direct circle kernel.
- 2026-03-13: F4 completed.
  - added `bench/julia/cuda/profile_interpolation.jl` to gate any future tiling/shared-memory work on actual profile evidence.
  - current candidate set remains interpolation-family kernels only.
- 2026-03-13: follow-up optimization slice implemented.
  - host-side CUDA interpolation profiling now synchronizes once after the profiled loop, so the host profile reflects launch/wrapper overhead instead of per-iteration blocking.
  - common full-frame `prop_end!` now uses direct quadrant copies / broadcasts for both complex and intensity outputs.
  - circular aperture now has a centered-circle KA specialization in addition to the general shifted-circle path.
- 2026-03-13: CUDA steady-state benchmark reconciliation completed.
  - the benchmark driver now uses the standalone FP64 workload run as the single source of truth for the compatibility-facing CUDA steady-state row.
  - the standard CUDA steady-state report is now emitted as an alias of the FP64 workload report rather than a second independent run.
  - summary reporting now warns if the standard CUDA steady-state and standalone FP64 reports diverge materially, which protects against stale manual runs.
- 2026-03-13: GPU follow-up slice implemented after benchmark reconciliation.
  - `prop_end!` full-frame dispatch was temporarily changed to a CUDA-specific single-kernel shifted-copy path, but that regressed on the RTX 3050 Ti and was removed.
  - centered standard circular aperture was temporarily changed to a box-specialized centered-aperture path, but that also regressed on the RTX 3050 Ti and was removed.
  - the retained change from this investigation is the benchmark conclusion itself: for these kernels on this hardware, fewer launches was not worth the worse memory access pattern in `prop_end!`, and multi-kernel box decomposition was not worth the added launch overhead for centered circular aperture.
- 2026-03-13: CUDA wavefront microbenchmark methodology unified.
  - the shared helper in `bench/common/cuda_wavefront_kernel_cases.jl` now defines the common setup for `prop_qphase`, `prop_ptp`, `prop_wts`, `prop_stw`, `prop_circular_aperture`, and `prop_end_mutating`.
  - both the supported-kernel report and the precision-split report now benchmark those kernels through the same state-restore and warmup path.
  - remaining discrepancies between those reports should now be interpreted as measurement variance or true kernel differences, not benchmark scaffolding drift.
