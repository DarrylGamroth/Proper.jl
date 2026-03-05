# Compatibility Decisions

This log records decisions when Python 3.3.4, MATLAB 3.3.1, and manual intent diverge.

## D-0001: Reference Implementation and Parity Baseline
- Date: 2026-03-04
- Status: Accepted
- Context:
  - MATLAB 3.3.1 appears closer to intent in some routines, but is not executable in this environment.
  - Python 3.3.4 is executable and can be used for automated numeric comparisons.
- Decision:
  - Use Python 3.3.4 as the primary executable baseline for parity/regression tests.
  - Use MATLAB 3.3.1 and the PROPER manual as semantic references to identify likely translation defects and intended behavior.
  - Implement dual compatibility modes in Julia:
    - `compat_mode = :python334` for strict Python parity.
    - `compat_mode = :corrected` for behavior corrected using MATLAB/manual-backed rationale.
- Consequences:
  - Automated golden-data generation and CI parity checks run against Python outputs.
  - Divergences from Python must be documented in this file with justification and tests.
  - Any corrected behavior must remain opt-in unless explicitly promoted by project decision.

## D-0002: FITS Library Selection
- Date: 2026-03-04
- Status: Accepted
- Context:
  - The port requires robust FITS read/write support for maps, headers, and parity artifacts.
  - A concrete library choice avoids interface churn during implementation.
- Decision:
  - Use `FITSIO.jl` (https://github.com/JuliaAstro/FITSIO.jl) for FITS file handling in Julia.
- Consequences:
  - `FITSIO.jl` is treated as a required dependency for FITS-related modules (`prop_fits_read`, `prop_fits_write`, `prop_writemap`, `prop_readmap`, error-map workflows).
  - Any alternate FITS backend would require a new compatibility decision and migration plan.

## D-0003: Default Plotting Package
- Date: 2026-03-04
- Status: Accepted
- Context:
  - Ported examples and demos need a consistent default plotting interface.
  - A single default avoids fragmented plotting code across examples.
- Decision:
  - Use `Plots.jl` as the default plotting package for PROPER Julia examples/demos.
- Consequences:
  - Example ports should target `Plots.jl` APIs by default.
  - Alternate plotting stacks (e.g., Makie/PyPlot) are optional and out of default scope unless explicitly requested.

## D-0004: Package and Module Naming
- Date: 2026-03-04
- Status: Accepted
- Decision:
  - Canonical package/module name is `Proper` (idiomatic Julia).
  - Preserve PROPER familiarity by keeping `prop_*` function names unchanged.
  - Provide a lightweight compatibility alias path for `proper` naming where practical.

## D-0005: Public API Compatibility Surface
- Date: 2026-03-04
- Status: Accepted
- Decision:
  - Keep public routine names as `prop_*` for familiarity.
  - Accept uppercase compatibility keywords and idiomatic lowercase keyword aliases.
  - Preserve Python return conventions for key entry points (`prop_run`, `prop_end`, `prop_run_multi`), while using concrete Julia types.

## D-0006: Default Compatibility Mode
- Date: 2026-03-04
- Status: Accepted
- Decision:
  - Default `compat_mode = :python334`.
  - Provide `compat_mode = :corrected` for MATLAB/manual-backed fixes.
  - Promote corrected behavior to default only through explicit future decision.

## D-0007: Numerical Convention Contract
- Date: 2026-03-04
- Status: Accepted
- Decision:
  - Match Python 3.3.4 normalization and centering behavior in `:python334` mode.
  - Freeze coordinate/pixel-center conventions per PROPER behavior and test them explicitly (even/odd sizes).
  - Expose all numerically meaningful deviations only via explicit compat-mode switches.

## D-0008: GPU Scope for Phase 1
- Date: 2026-03-04
- Status: Accepted
- Decision:
  - GPU-ready in Phase 1: core field algebra and FFT-based propagation kernels.
  - CPU-first in early phases: FITS I/O, complex interpolation/rotation/DM fitting paths.
  - Enforce `CUDA.allowscalar(false)` in GPU tests.

## D-0009: Backend Traits and Interfaces
- Date: 2026-03-04
- Status: Accepted
- Decision:
  - Define trait-dispatched backends for FFT, interpolation/resampling, and RNG.
  - Keep map I/O on `FITSIO.jl` with explicit host/device adaptation boundaries.
  - Require stable internal interfaces before heavy kernel porting begins.

## D-0010: Golden Data Workflow
- Date: 2026-03-04
- Status: Accepted
- Decision:
  - Generate baselines from Python 3.3.4 using reproducible scripts.
  - Store compact baseline artifacts + provenance metadata under `test/parity/`.
  - Keep full-size optional artifacts out of the default git path unless explicitly approved.

## D-0011: RNG Determinism Policy
- Date: 2026-03-04
- Status: Accepted
- Decision:
  - Library runtime default remains non-deterministic (familiar behavior).
  - Tests/parity runs must use explicit seeds.
  - For cross-backend parity, generate random phases deterministically on CPU and transfer when needed.

## D-0012: `prop_run_multi` Execution Semantics
- Date: 2026-03-04
- Status: Accepted
- Decision:
  - Default execution model: shared-memory threading (Julia `Threads`).
  - Optional distributed execution may be added behind explicit keyword.
  - Output ordering must match input ordering deterministically.

## D-0013: Performance Targets
- Date: 2026-03-04
- Status: Accepted
- Decision:
  - Core elementwise operations (`prop_add_phase`, `prop_multiply`, `prop_divide`) should be allocation-free in steady state.
  - Propagation kernels should use preallocated workspace in benchmark paths.
  - Initial acceptance target: no worse than 1.5x Python baseline on CPU for representative workloads, then optimize downward.

## D-0014: CI and Support Matrix
- Date: 2026-03-04
- Status: Accepted
- Decision:
  - Minimum supported Julia: 1.10.
  - Required CI: Linux CPU (parity + unit).
  - Optional CI: additional Julia version, macOS smoke, nightly GPU parity job.

## D-0015: Legal and Fixture Policy
- Date: 2026-03-04
- Status: Accepted
- Decision:
  - Keep source licensing files and attribution intact.
  - Commit only fixtures/artifacts that are legally redistributable and size-appropriate.
  - Keep large/generated parity products outside default repository history unless explicitly approved.

## D-0016: Scope for MATLAB-Only Utilities
- Date: 2026-03-04
- Status: Accepted
- Decision:
  - Implement `prop_hex_aperture` and `prop_version` for user familiarity.
  - Treat MATLAB internals (`propcommon`, `prop_free_threads`) as non-portable implementation details.
  - Treat plotting helpers (`figarray`, `figplace`) as out-of-core optional examples.

## D-0017: `compat_mode` Evaluation Policy
- Date: 2026-03-04
- Status: Accepted
- Context:
  - Repeated `compat_mode` branching inside inner loops adds overhead and increases behavioral complexity.
  - Mode decisions should be deterministic and easy to audit.
- Decision:
  - Accept `compat_mode` only as a keyword in context/config constructor(s).
  - Resolve `compat_mode` exactly once during context/config construction.
  - Convert mode to an internal typed policy and pass that policy through call chains.
  - Do not branch on mode inside hot-path inner loops.
- Consequences:
  - Public `prop_*` APIs operate on pre-built context/config objects and do not re-accept `compat_mode`.
  - Kernel methods can dispatch on policy type, reducing runtime branching overhead.
  - Raw mode symbols are constructor-only inputs and are not propagated through runtime call chains.
  - Removing `compat_mode` later is localized to constructor/context and policy mapping code.

## D-0018: API Layering With Mutating Kernels
- Date: 2026-03-04
- Status: Accepted
- Decision:
  - Keep user-facing `prop_*` routines for familiarity.
  - Implement performance-critical internals as mutating `*_!` kernels.
  - Use `prop_*` wrappers for argument normalization and contract checks.

## D-0019: Typed Run Context
- Date: 2026-03-04
- Status: Accepted
- Decision:
  - Introduce a typed run context carrying compat policy, backend traits, RNG, and workspace handles.
  - Pass context explicitly through call chains instead of using module-global state.

## D-0020: One-Time Option Normalization
- Date: 2026-03-04
- Status: Accepted
- Decision:
  - Parse/normalize compatibility keywords and PASSVALUE inputs once at API boundary.
  - Convert dynamic containers to typed option structs or `NamedTuple`s before hot-path execution.

## D-0021: Preallocated Workspace Policy
- Date: 2026-03-04
- Status: Accepted
- Decision:
  - Use reusable workspace buffers for propagation and interpolation routines.
  - Avoid per-call temporary allocations in steady-state workloads.

## D-0022: FFT Plan Caching Policy
- Date: 2026-03-04
- Status: Accepted
- Decision:
  - Use cached FFT plans (CPU/GPU) via backend traits.
  - Do not build FFT plans inside inner loops.

## D-0023: Type-Stable Core Structs
- Date: 2026-03-04
- Status: Accepted
- Context:
  - High-frequency propagation kernels are sensitive to inference loss.
  - Type instability causes allocations and hidden dynamic dispatch.
- Decision:
  - Core structs use concrete/parametric fields only.
  - Avoid abstractly typed fields in hot-path state containers.
  - Avoid `Any`-typed intermediates in hot-path kernels.
- Consequences:
  - New hot-path code must include type-stability validation in tests/bench checks.
  - Abstract field usage requires explicit justification and should remain out of inner loops.

## D-0024: Optional Input Representation
- Date: 2026-03-04
- Status: Accepted
- Decision:
  - Accept legacy dynamic input forms for compatibility.
  - Normalize to typed options / `NamedTuple` once, then dispatch on typed values.

## D-0025: RNG Injection Policy
- Date: 2026-03-04
- Status: Accepted
- Decision:
  - Support explicit RNG injection through context/options.
  - Keep runtime default non-deterministic, but require seeded RNG for tests/parity.

## D-0026: Backend-Agnostic First, Specialize By Profiling
- Date: 2026-03-04
- Status: Accepted
- Decision:
  - Implement backend-agnostic `AbstractArray` paths first.
  - Add backend-specialized methods only where profiling shows measurable need.

## D-0027: Internal Namespace Organization
- Date: 2026-03-04
- Status: Accepted
- Decision:
  - Keep one-to-one `src/prop_*.jl` compatibility files.
  - Place reusable internals in structured submodules (e.g., `src/core/`) and call them from compatibility files.

## D-0028: Dynamic Dispatch Minimization
- Date: 2026-03-04
- Status: Accepted
- Context:
  - Repeated runtime dispatch in tight loops harms throughput and predictability.
  - Dispatch is still useful at API boundaries and for coarse-grained backend selection.
- Decision:
  - Limit dynamic dispatch to API/adapter boundaries.
  - Inner-loop kernels should run on concretely inferred types and policy/trait-dispatched methods.
  - Mode/backend/options branching should be resolved before entering hot loops.
- Consequences:
  - Performance checks should include dispatch/inference inspection for hot kernels.
  - New features should prefer pre-resolved typed contexts over runtime symbol/dict branching.

## D-0029: Benchmark Methodology And TTFx Policy
- Date: 2026-03-04
- Status: Accepted
- Context:
  - We need fair Python-vs-Julia runtime comparisons and actionable optimization data.
  - Julia compilation latency (TTFx) can dominate first-run timings and obscure steady-state performance.
- Decision:
  - Primary performance comparisons use steady-state timings and explicitly exclude Julia TTFx.
  - Julia benchmarks must run after warmup/compilation and in a persistent process.
  - Track TTFx separately as an independent metric for startup/precompile work.
  - Use benchmark results to drive precompile statements/workloads.
  - Do not require `PrecompileTools.jl` by default; rely on modern Julia startup behavior unless TTFx data indicates targeted precompile directives are needed.
- Consequences:
  - Benchmark reports include two sections:
    - steady-state throughput/latency (Python vs Julia)
    - Julia cold-start/TTFx diagnostics (not mixed into steady-state comparison)
  - Precompile improvements are evaluated against the separate TTFx metric.

## D-0030: Final Parity Goal Definition
- Date: 2026-03-04
- Status: Superseded-in-part by D-0031
- Context:
  - Exact bit-level Python equivalence is not required when differences are numerically negligible in deep-null regions.
  - MATLAB 3.3.1 is a useful style/semantic reference even though runtime comparison is unavailable in this environment.
- Decision:
  - Final target is:
    - algorithm structure semantically comparable to upstream PROPER behavior,
    - numerical behavior physics-equivalent to Python 3.3.4 for executable parity workflows.
  - Parity acceptance must use both relative and absolute error criteria (not relative-only) for high-contrast/null outputs.
- Consequences:
  - Remaining parity closure focuses on physically meaningful agreement and documented threshold criteria.
  - MATLAB remains a semantic reference; Python remains the executable numeric baseline.
  - Implementation style guidance is defined by D-0031 (idiomatic Julia, not MATLAB style mirroring).

## D-0031: Implementation Style Clarification
- Date: 2026-03-04
- Status: Accepted
- Context:
  - The project should not target MATLAB coding style directly.
  - We still want algorithmic comparability with upstream implementations while preserving Julia quality.
- Decision:
  - Prefer idiomatic Julia implementation patterns (multiple dispatch, concrete types, traits, type-stable kernels).
  - Validate that algorithms are semantically similar to upstream PROPER implementations.
  - Keep parity target as physics-equivalent behavior to Python 3.3.4 (with documented absolute/relative acceptance criteria).
- Consequences:
  - MATLAB is treated as semantic reference material only, not a style template.
  - Code review criteria prioritize Julia idioms and performance constraints over direct stylistic mirroring.

## D-0032: Current Corrected-Mode Divergence Policy
- Date: 2026-03-04
- Status: Accepted
- Context:
  - Phase 8 exit criteria require any remaining `:corrected` deltas to be intentional and documented.
  - At present, no active algorithmic behavior divergence is implemented between `:python334` and `:corrected`.
- Decision:
  - Treat `:corrected` as behavior-equivalent to `:python334` until specific corrected algorithms are intentionally introduced.
  - Any future corrected-mode divergence must be:
    - explicitly documented in `docs/compat_decisions.md`,
    - covered by tests for both modes,
    - reflected in parity threshold rationale where relevant.
- Consequences:
  - Current phase-8 parity closure is evaluated against `:python334` thresholds with no undocumented corrected-mode deltas.
  - Corrected-mode behavior changes remain controlled and auditable.
