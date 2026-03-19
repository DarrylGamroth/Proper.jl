# PROPER Python 3.3.4 -> Julia Port Plan

## Goals
- Port `../proper_v3.3.4_python` to Julia with an idiomatic API.
- Keep a one-to-one filename mapping between Python and Julia sources.
- Use `AbstractArray`-based implementations so CPU arrays and GPU arrays are both supported.
- Port all upstream examples and validate Julia vs Python numerical agreement.
- Keep filename parity for traceability, but do not do line-by-line transliteration; use Julia-native types, dispatch, and trait-based backend selection.

## Porting Principles (Julia-First Internals)
- Filename mapping is an organization constraint, not an implementation constraint.
- Preserve externally visible behavior and naming where needed for compatibility.
- Refactor internals into typed kernels and multiple-dispatch methods.
- Prefer trait-based backend routing (`Array`, `CuArray`, future backends) over runtime flag toggles.
- Prefer composable kernels (`KernelAbstractions.jl`, `AcceleratedKernels.jl`) over Python-style global switches and ad hoc branching.

## Package Review Summary
- Source size: 87 Python modules in `proper/`, 23 Python example scripts in `proper/examples/`, 3 C kernels, 1 FITS data asset.
- Core architecture today:
  - `proper/__init__.py` re-exports nearly all routines and stores mutable global run state.
  - `WaveFront` object (`prop_wavefront.py`) holds optical state and array payload.
  - `prop_run`/`prop_run_multi` dynamically import prescriptions and rely on mutable module globals.
  - FFT backend is selected via side effects (`numpy`, `pyfftw`, Intel MKL wrappers).
  - Interpolation and map transforms use SciPy or optional C shared libraries.
- Main migration pressure points for Julia:
  - Global mutable state and keyword-string switches are not thread-safe and not Julia-idiomatic.
  - Many routines assume `numpy.ndarray` and `Float64/Complex128`; these need to be generalized to `AbstractArray`.
  - SciPy-specific interpolation/rotation paths need Julia-native CPU/GPU-compatible kernels.
  - Several routines use dynamic `exec/eval`, which should be replaced by explicit generated algebra or dispatch.

## Critical Findings In Upstream Python (Port Should Intentionally Address)
1. `prop_psd_errormap.py` calls `proper.prop_use_fftw()` inside compute path and checks `== 1`, but `prop_use_fftw` returns `None`; FFTW branch is effectively unreachable and function has unwanted side effects. (`proper/prop_psd_errormap.py:185`)
2. `prop_errormap.py` uses `proper_switch_set` (undefined symbol), causing runtime `NameError` in some keyword combinations. (`proper/prop_errormap.py:92`)
3. Historical Python 3.3.4 issue: `prop_resamplemap.py` used `xshift` for both axes. The local parity baseline now patches this to use `yshift` on the Y axis, matching MATLAB intent. (`proper/prop_resamplemap.py:67`)
4. `prop_end.py` uses float slice indices for `EXTRACT` under Python 3 (`ny/2`, `nx/2`), which is invalid for slicing. (`proper/prop_end.py:52`)
5. `prop_state.py` assigns loaded state to local `wf` only; caller's wavefront object is not updated. (`proper/prop_state.py:49`)

## MATLAB v3.3.1 Cross-Review (Secondary Baseline)
- Reviewed `../proper_v3.3.1_matlab` as an additional reference for numerical intent.
- MATLAB-only core files (not in Python 3.3.4): `propcommon.m`, `prop_free_threads.m`, `prop_version.m`, `prop_readfits.m`, `prop_hex_aperture.m`.
- Python-only core files (not in MATLAB 3.3.1): package/export wrappers and backend-specific FFT/interpolation integration (`prop_fftw*`, `prop_ffti*`, `prop_use_fft*`, `prop_compile_c`, `prop_cubic_conv`, `prop_wavefront`, etc.).
- Example differences: MATLAB includes plotting helpers `figarray.m`, `figplace.m`; Python includes package `__init__.py`.
- Cross-check outcomes on drift-prone routines:
  - MATLAB `prop_resamplemap` correctly applies `yshift` separately. (`proper_v3.3.1_matlab/prop_resamplemap.m:49`)
  - MATLAB `prop_end` uses explicit integer center slicing for `extract`. (`proper_v3.3.1_matlab/prop_end.m:60`)
  - MATLAB `prop_state` returns an updated beam struct, avoiding Python's local-rebind issue. (`proper_v3.3.1_matlab/prop_state.m:54`)
  - MATLAB `prop_psd_errormap` computes directly via FFT path without calling backend-toggle setup routines in-band. (`proper_v3.3.1_matlab/prop_psd_errormap.m:222`)
- Note: MATLAB and Python both appear to share an in-place rotation-variable reuse in `prop_psd_errormap` coordinate setup; treat as a compatibility quirk pending manual/IDL confirmation. (`proper_v3.3.1_matlab/prop_psd_errormap.m:196`)

## Julia Target Architecture (Idiomatic + GPU-Capable)
### 1. Public API and Module Layout
- Keep one Julia file per Python source file, same basename (except package entrypoint).
- Top-level module file (`src/Proper.jl`) includes and re-exports per-file functions.
- Preserve PROPER function names (`prop_begin`, `prop_lens`, etc.) for compatibility.
- Add optional higher-level idiomatic aliases later (non-breaking).
- Use compatibility wrappers at file boundaries, with shared typed internals (not duplicated Python-structure logic in every file).

### 2. State Model
- Replace module-global mutable state with explicit structs:
  - `ProperConfig` for runtime options (verbosity, FFT backend preference, antialias settings).
  - `ProperRuntime` for per-run mutable bookkeeping (table arrays, save-state metadata).
  - `WaveFront{T,A<:AbstractMatrix{Complex{T}}}` for optical field and propagation state.
- Thread/process safety comes from passing context/state explicitly.

### 3. Dispatch and Trait Strategy
- Introduce backend traits (e.g., `BackendStyle`) selected from array/type inputs rather than mutable globals.
- Dispatch core operations by trait:
  - FFT operations (`fft2!`, `ifft2!`)
  - interpolation/resampling operations
  - elementwise phase/amplitude transforms
- Keep high-level `prop_*` APIs stable while routing to specialized kernels via dispatch.

### 4. AbstractArray and GPU Strategy
- All array arguments should be typed as `AbstractArray` / `AbstractMatrix` and preserve input backend via `similar`/`eltype`.
- No hardcoded `Array` allocations in hot paths; allocate through helper constructors that preserve array backend.
- FFT dispatch via backend traits:
  - CPU default: `FFTW.jl` (through `AbstractFFTs`).
  - GPU path: `CUDA.jl` (`CUDA.CUFFT`) when arrays are `CuArray`.
  - Optional MKL path if needed.
- Interpolation/rotation/resampling:
  - Phase 1: robust CPU implementation in pure Julia.
  - Phase 2: GPU-capable kernels using `KernelAbstractions.jl` and `AcceleratedKernels.jl` for `prop_cubic_conv`, `prop_szoom`, map resampling, and rotation.

### 5. Numerical and Performance Guidance
- Keep computations inside concrete, type-stable functions (no untyped globals).
- Prefer in-place `!` internal kernels for hot paths to reduce allocations.
- Avoid dynamic `eval/exec`; replace with explicit formula evaluation and precomputed terms.
- Add `@inbounds`/`@views` only after correctness parity is established.

### 6. Compatibility Policy
- Default behavior should match Python outputs for existing examples.
- For known Python bugs above, either:
  - patch the executable baseline and regenerate artifacts, or
  - document and test any intentional Julia-side divergence.

## Phase 0: Preflight Decisions (Required Before Phase 1)
Purpose: freeze compatibility, architecture, and validation contracts so file-by-file porting does not thrash.

### 0.1 API and Compatibility Contract
- Decide preserved public surface (`prop_*` function names, keyword names/casing, return-value shapes/types).
- Decide which known Python bugs are baseline-patched versus intentionally retained.
- Decide whether keyword parsing will support both uppercase compatibility style and idiomatic Julia keyword aliases.

### 0.2 Naming and Packaging
- Confirm module/package naming (`Proper` canonical with compatibility alias, or strict `proper`).
- Confirm long-term directory policy while keeping one-to-one filename mapping.

### 0.3 Numerical Contract
- Freeze FFT normalization conventions and forward/inverse scaling.
- Freeze centering conventions (`prop_shift_center` / `fftshift` semantics).
- Freeze coordinate and pixel-center conventions used by aperture/mask/interpolation routines.
- Freeze unit conventions (meters, microns, radians, arcsec) and conversion locations.

### 0.4 Backend Trait Contract
- Define trait interfaces and required methods for:
  - FFT backend
  - interpolation/resampling backend
  - random number generation backend
  - map IO backend
- Define fallback order (pure Julia CPU baseline, then optional accelerated paths).

### 0.5 GPU Scope for Phase 1
- Explicitly mark which routines must be GPU-ready in initial phases vs CPU-only with planned follow-up.
- Enforce no scalar GPU indexing policy in tests.

### 0.6 Determinism and Reproducibility
- Decide deterministic policy for randomized routines (`prop_psd_errormap` and any random phase generation).
- Define seed strategy and artifact metadata requirements for parity fixtures.

### 0.7 Regression Harness Contract
- Define golden-data workflow (Python baseline generation, storage format, provenance metadata).
- Define parity metrics and thresholds by output class (intensity vs complex field, CPU vs GPU).
- Define acceptable tolerances for each example family.

### 0.8 Performance Contract
- Set baseline allocation and runtime budgets for core kernels:
  - `prop_ptp`, `prop_wts`, `prop_stw`, `prop_qphase`, `prop_lens`, `prop_propagate`.
- Define benchmark environments and pass/fail criteria.

### 0.9 Parallel Execution Contract
- Decide semantics for Julia equivalent of `prop_run_multi` (`Threads` vs `Distributed`).
- Decide determinism expectations for parallel execution and random workloads.

### 0.10 Dependency, CI, and Release Constraints
- Freeze minimum Julia version and supported CUDA/toolchain versions.
- Classify dependencies as required vs optional (`KernelAbstractions`, `AcceleratedKernels`, `CUDA`, `MKL`).
- Define CI matrix (quick CPU smoke, full CPU parity, optional GPU parity).
- FITS handling dependency is fixed to `FITSIO.jl` (see `D-0002` in `docs/compat_decisions.md`).
- Default plotting dependency for examples/demos is `Plots.jl` (see `D-0003` in `docs/compat_decisions.md`).

### 0.11 Legal and Data Governance
- Confirm `LegalStuff.txt` constraints for redistribution of bundled assets and generated reference artifacts.
- Define what fixtures/results may be checked into version control.

### 0.12 Cross-Implementation Precedence Rule
- Define authoritative source when implementations disagree:
  - Python 3.3.4 for public API and file-level compatibility.
  - MATLAB 3.3.1 for numerical intent when Python behavior is clearly a translation defect.
  - PROPER manual as tie-breaker when both codebases are ambiguous.
- Track each divergence decision in a dedicated compatibility log (recommended: `docs/compat_decisions.md`).
- Accepted decisions are tracked in `docs/compat_decisions.md` (`D-0001` through `D-0029`).

Exit criterion for Phase 0: each subsection above is marked “decided” in project tracking and reflected in docs/tests; accepted decisions are recorded in `docs/compat_decisions.md` (`D-0001` through `D-0029`).

## Implementation Phases
### Phase 0: Preflight Contracts (Required)
- Goals:
  - Keep decisions in `docs/compat_decisions.md` current and enforceable (`D-0001` through `D-0029` currently accepted).
  - Freeze API, numerics, trait interfaces, and test/provenance rules.
- Work items:
  - Maintain compatibility decision log for any new divergence/design change.
  - Write short contract docs:
    - `docs/api_contract.md`
    - `docs/numerics_contract.md`
    - `docs/backend_traits.md`
    - `docs/parity_harness_contract.md`
  - Define required CI jobs and minimum Julia/toolchain versions.
- Deliverables:
  - All phase-0 decisions marked `Accepted` or `Rejected`.
  - Contract docs present and cross-linked from `PORTING_PLAN.md`.
- Exit criteria:
  - No unresolved phase-0 decision blocks implementation.
  - A contributor can implement any module without re-litigating contracts.

### Phase 1: Foundation And Infrastructure
- Goals:
  - Create stable project scaffolding and runtime model.
  - Establish compatibility-layer patterns used by all subsequent ports.
- Work items:
  - Build full source file skeleton with one-to-one filename mapping and deterministic `include(...)` order.
  - Add internal namespace layout for reusable core internals (e.g., `src/core/`) while keeping compatibility files in `src/prop_*.jl`.
  - Implement core types:
    - `WaveFront{T,A<:AbstractMatrix{Complex{T}}}`
    - runtime/config structs for formerly global state
    - typed run context carrying backend/RNG/workspace
  - Implement keyword normalization (compat uppercase + idiomatic aliases).
  - Normalize optional inputs once at API boundary into typed options.
  - Implement backend trait stubs and default CPU backend wiring.
  - Define mutating internal kernel naming/pattern (`*_!`) and wrapper policy (`prop_*`).
  - Define workspace and FFT-plan caching interfaces.
  - Add test scaffolding:
    - unit test layout by module family
    - parity harness entrypoints (Python invocation + metadata capture)
- Deliverables:
  - Compile-able package skeleton with stable module boundaries.
  - Baseline tests for constructors/getters, keyword parsing, and trait dispatch.
- Exit criteria:
  - Package loads cleanly.
  - Core structs and trait APIs are stable enough to support feature phases.

### Phase 2: Core Propagation Kernel
- Goals:
  - Port end-to-end wavefront initialization, propagation, and termination.
  - Establish performance patterns for hot paths.
- Work items:
  - Port and verify:
    - `prop_begin`, `prop_end`
    - `prop_propagate`, `prop_select_propagator`
    - `prop_ptp`, `prop_stw`, `prop_wts`, `prop_qphase`
    - `prop_lens`
  - Implement shared internal helpers for centered frequency grids and reusable workspace.
  - Add behavior checks against the patched Python baseline where needed.
  - Benchmark core kernels for allocations and throughput.
  - Add type-inference/dispatch checks for hot kernels.
- Deliverables:
  - Core propagation path operational on CPU with deterministic tests.
  - Benchmark baselines recorded for hot kernels.
- Exit criteria:
  - `simple_prescription` and `simple_telescope` run with parity thresholds.
  - Core kernel allocation targets are met in steady-state benchmark paths.

### Phase 3: Apertures, Obscurations, And Field Algebra
- Goals:
  - Port all geometric aperture/obscuration primitives and map application ops.
- Work items:
  - Port geometry/mask builders:
    - `prop_ellipse`, `prop_rectangle`, `prop_polygon`, `prop_irregular_polygon`, `prop_rounded_rectangle`
    - circular/elliptical/rectangular aperture and obscuration wrappers
  - Port field algebra operations:
    - `prop_multiply`, `prop_divide`, `prop_add_phase`, `prop_add_wavefront`
  - Validate anti-aliasing behavior and keyword switches (`NORM`, `DARK`, rotations).
- Deliverables:
  - Complete mask/algebra module family with compatibility tests.
- Exit criteria:
  - Shape/mask primitives pass pixel-level regression tests on representative grids.
  - Field algebra routines are type-stable and allocation-aware.

### Phase 4: FITS, Map Pipeline, And Resampling/Interpolation
- Goals:
  - Port all map IO and transforms used by PSD/error-map and optics workflows.
- Work items:
  - Port FITS modules with `FITSIO.jl`:
    - `prop_fits_read`, `prop_fits_write`, `prop_writemap`, `prop_readmap`
  - Port map/error generation:
    - `prop_errormap`, `prop_psd_errormap`, `prop_resamplemap`
  - Port interpolation/rotation stack:
    - `prop_magnify`, `prop_szoom`, `prop_cubic_conv`, `prop_rotate`
  - Implement CPU baseline kernels first, then accelerated kernels via `AcceleratedKernels.jl` / `KernelAbstractions.jl`.
  - Add deterministic RNG hooks for parity tests in randomized map generation.
- Deliverables:
  - Full map workflow available and tested against patched Python baseline.
- Exit criteria:
  - `psdtest` and map round-trip tests pass parity thresholds.
  - FITS headers/units semantics match contract.

### Phase 5: Zernikes And Segmented Optics
- Goals:
  - Port aberration fitting/generation stack and segmented-aperture functionality.
- Work items:
  - Port:
    - `prop_noll_zernikes`, `prop_zernikes`, `prop_fit_zernikes`, `prop_print_zernikes`
    - `prop_hex_zernikes`, `prop_hex_wavefront`
  - Replace Python dynamic `exec/eval` behavior with explicit generated/precomputed terms.
  - Add numerical conditioning checks for fitting routines.
- Deliverables:
  - Zernike and segmented optics modules with robust tests and docs.
- Exit criteria:
  - Zernike generation/fitting parity tests pass on unobscured and obscured cases.
  - Segmented-aperture examples produce expected fields/maps.

### Phase 6: DM, Multi-Run, And Save-State
- Goals:
  - Port deformable mirror path, parallel prescription execution, and state lifecycle.
- Work items:
  - Port:
    - `prop_dm`, `prop_fit_dm`
    - `prop_run_multi` and supporting execution helpers
    - state modules (`prop_init_savestate`, `prop_state`, `prop_savestate`, `prop_end_savestate`, `prop_is_statesaved`)
  - Implement deterministic ordering and seeded behavior in parallel runs.
  - Resolve state-file format policy and document compatibility guarantees.
- Deliverables:
  - Multi-run and DM workflows integrated into parity harness.
- Exit criteria:
  - `multi_example`, `testmulti1`, `testmulti2`, and DM-driven examples pass parity checks.
  - Save/restore behavior is validated in integration tests.

### Phase 7: Examples, Parity Harness, And Release Readiness
- Goals:
  - Complete user-facing example parity and release-quality validation.
- Work items:
  - Port all 23 Python examples one-to-one (plus accepted MATLAB-only additions if approved).
  - Finalize parity harness:
    - Python baseline generation scripts
    - artifact provenance capture
    - tolerance reports per backend/mode
  - Finalize benchmarking harness:
    - steady-state Python-vs-Julia benchmark suite
    - separate Julia TTFx/cold-start benchmark suite
    - benchmark artifact export (JSON/CSV + markdown summaries)
  - Finalize CI matrix:
    - required CPU jobs
    - optional GPU jobs
  - Write migration notes and compatibility guide for users.
- Deliverables:
  - Fully ported examples, parity reports, benchmark reports, and CI pipelines.
- Exit criteria:
  - Acceptance criteria section is fully satisfied.
  - Project is ready for tagged pre-release.

### Phase 8: Full Parity Closure (Python 3.3.4 Baseline)
- Goals:
  - Eliminate placeholder/fallback behaviors in physics-critical modules.
  - Achieve numerical parity against executable Python baseline across full example suite.
- Work items:
  - Replace fallback implementations with parity implementations for:
    - propagation internals (`prop_propagate`, `prop_ptp`, `prop_stw`, `prop_wts`, `prop_select_propagator`)
    - PSD/map synthesis and transforms (`prop_psd_errormap`, interpolation/cubic-conv paths)
    - Zernike/fit stack (`prop_zernikes`, `prop_noll_zernikes`, `prop_fit_zernikes`, related helpers)
    - polygon/segmented optics paths (`prop_polygon`, `prop_irregular_polygon`, `prop_hex_*`, rounded geometry)
  - Run parity harness on all 23 examples against the patched Python baseline.
  - Add per-module parity fixtures and failure triage reports with decision-log links.
  - Close parity gaps by updating implementation and/or explicitly documenting accepted divergence.
- Deliverables:
  - No physics-critical module relies on placeholder/fallback behavior.
  - Full parity report for all examples with provenance and threshold outcomes.
- Exit criteria:
  - All example parity thresholds are met against the patched Python baseline.
  - Any remaining intentional divergences are documented in `docs/compat_decisions.md`.
  - No unresolved high-severity parity gaps remain.

### Phase 9: MATLAB Semantic Reconciliation And Final Validation
- Goals:
  - Reconcile known Python translation defects using MATLAB/manual semantics where applicable.
  - Prepare final release-quality validation matrix and migration notes.
- Work items:
  - Audit known disagreement hotspots against MATLAB source/manual intent:
    - `prop_resamplemap` shift semantics
    - `prop_end` extract indexing semantics
    - `prop_state` restore semantics
    - `prop_psd_errormap` backend-toggle/side-effect behavior
  - Promote or retain corrected behavior per accepted decisions and compat-mode policy.
  - Add explicit MATLAB-semantic regression tests (non-executable reference checks).
  - Finalize migration guide for users moving from Python/MATLAB to Julia.
- Deliverables:
  - Final semantic reconciliation report (Python parity + MATLAB/manual rationale).
  - Release checklist signed off with decision references.
- Exit criteria:
  - Python parity baseline is green and MATLAB/manual-backed corrections are documented and tested.
  - Release notes clearly describe compatibility guarantees and corrected-mode behavior.

## Refactor Track (Post-Audit, Parity-Preserving)
Context: `docs/archive/JULIA_IMPLEMENTATION_AUDIT_2026-03-04.md` immediate checklist is complete. This track closes remaining architecture/performance findings without changing user-facing behavior.

### R1: Typed Option Boundaries (Remove Dynamic Keyword Hot Paths)
- Goals:
  - Eliminate repeated `kwargs`/`haskey` branching inside compute paths.
  - Normalize option parsing once at API boundary (`D-0020`, `D-0028`).
- Work items:
  - Add typed option structs + normalizers for:
    - propagation family (`prop_propagate`, `prop_select_propagator`, `prop_lens`)
    - map/error family (`prop_errormap`, `prop_psd_errormap`, `prop_rotate`, `prop_magnify`)
    - geometry family (`prop_ellipse`, `prop_rectangle`, `prop_polygon`, `prop_irregular_polygon`)
  - Keep compatibility keywords in wrappers only; kernels accept typed options.
- Exit criteria:
  - No runtime keyword dictionary logic in inner loops of designated hot kernels.
  - `@code_warntype` clean on normalized kernel entry points.

### R2: Trait-Driven Kernel Routing (CPU/GPU-Ready)
- Goals:
  - Move backend/policy selection to trait dispatch points instead of runtime branching.
  - Ensure backend preservation for `AbstractArray` inputs.
- Work items:
  - Define/extend internal dispatch entrypoints:
    - interpolation trait (`interp_backend(::A)` style) for resample/rotate/magnify/cubic-conv grid path
    - FFT trait usage in propagation kernels
  - Thread `RunContext` through internal kernels that currently infer backend implicitly.
  - Add explicit fallback hierarchy: generic CPU loop -> accelerated CPU/GPU kernel.
- Exit criteria:
  - Kernel backend selection visible as dispatch, not symbol/boolean branching.
  - GPU smoke tests execute for implemented trait-routed kernels with `CUDA.allowscalar(false)`.

### R3: Mutating Kernels + Workspace Reuse
- Goals:
  - Cut steady-state allocations in propagation/interpolation hot paths (`D-0018`, `D-0021`).
- Work items:
  - Introduce/complete `*_!` kernels for:
    - `prop_rotate!`, `prop_magnify!`, `prop_cubic_conv_grid!`
    - geometry mask fill kernels where full-frame temporaries are currently built
  - Add reusable workspace structs for repeated axis buffers, interpolation scratch, and FFT buffers.
  - Keep current `prop_*` non-mutating wrappers as compatibility adapters.
- Exit criteria:
  - Selected benchmark kernels meet defined allocation budgets in steady-state.
  - No regression in parity thresholds relative to current baseline.

### R4: WaveFront State Typing And Dispatch Simplification
- Goals:
  - Remove string/symbol state machine overhead in propagation control paths.
- Work items:
  - Replace `Symbol` propagation-state fields with compact typed selectors (`@enum` or equivalent small type tags).
  - Dispatch transition math on typed state selectors.
  - Preserve external behavior and table logging values.
- Exit criteria:
  - Propagation state transitions no longer rely on runtime symbol/string comparisons.
  - Transition kernels are inference-stable and parity tests remain green.

### R5: Performance/Inference Gates And Benchmark Expansion
- Goals:
  - Convert refactor expectations into enforced tests and benchmark gates (`D-0029`).
- Work items:
  - Expand inference/allocation tests beyond current gate set:
    - propagation core kernels
    - PSD/map synthesis hotspots
    - geometry/mask hotspots
  - Add benchmark matrix entries for refactored kernels and end-to-end examples.
  - Continue strict separation of steady-state vs cold-start/TTFx reporting.
- Exit criteria:
  - CI enforces inference/allocation checks on designated kernels.
  - Benchmark reports show no major regressions and document speedup deltas from refactor work.

### Refactor Completion Criteria
- Public `prop_*` API remains parity-stable.
- Hot paths are type-stable, dispatch-predictable, and allocation-bounded.
- Backend routing is trait-based and ready for progressive KA/AK rollout.
- All changes are covered by parity tests and documented decisions.

### KernelAbstractions/AcceleratedKernels Placement Plan
Goal: use `KernelAbstractions.jl` (KA) for custom index-heavy kernels and `AcceleratedKernels.jl` (AK) for high-level fused elementwise/reduction operations where coverage is sufficient.

| Area | Primary package | Candidate routines | Why |
| --- | --- | --- | --- |
| Grid-sampling interpolation | `KernelAbstractions.jl` | `prop_resamplemap!`, `prop_rotate!`, `prop_magnify!`, `prop_cubic_conv` grid path | gather/scatter indexing and custom boundary handling are clearer as explicit kernels |
| Geometry mask fill | `KernelAbstractions.jl` | `prop_ellipse`, `prop_rectangle`, `prop_irregular_polygon` fill loops | per-pixel tests with branchy geometry logic map naturally to custom kernels |
| PSD spectral-map synthesis | `KernelAbstractions.jl` + `AcceleratedKernels.jl` | `prop_psd_errormap` frequency-grid assembly + postprocessing | KA for coordinate/index transforms, AK for fused elementwise scaling/normalization |
| Field algebra and simple map transforms | `AcceleratedKernels.jl` | `prop_multiply`, `prop_divide`, `prop_add_phase`, selected map postprocess steps | mostly elementwise operations where high-level fused kernels can reduce temporaries |
| Reductions/statistics in map normalization | `AcceleratedKernels.jl` (if supported) else CPU fallback | RMS/peak normalization in PSD/error-map flows | use AK reductions where numerically consistent; fallback remains deterministic CPU path |

Integration policy:
- Keep generic CPU loops as correctness baseline.
- Route backend choice through traits/dispatch first, then select KA/AK implementation.
- Add KA/AK methods only when parity tests and steady-state benchmarks show clear wins.

## Example Parity Config Matrix (Coverage Plan)
| Area | Axis | Permutations | Evidence | Status | Notes |
| --- | --- | --- | --- | --- | --- |
| Example parity | Example script | 23/23 upstream examples | `test/example_parity/*.jl` + generated diff metrics | Gap | Must cover every example below |
| Backend | Array backend | CPU `Array`, GPU `CuArray` | Separate CI jobs + tolerance report | Gap | GPU first on core subset, then expand |
| FFT backend | FFT provider | FFTW default, optional MKL/CUFFT | backend-tagged regression report | Gap | Must keep same normalization conventions |
| Output mode | `prop_end` mode | intensity (`NOABS=false`), complex field (`NOABS=true`) | per-example assertions | Gap | Include talbot/talbot_correct |
| Multi-run | execution mode | single `prop_run`, parallel `prop_run_multi` | shape + value parity checks | Gap | Include PASSVALUE arrays |
| Error maps | map pipeline | read/resample/rotate/magnify/PSD | golden FITS + statistics checks | Gap | Include `psdtest`, coronagraph/telescope maps |
| State system | save/restore | save new state, load existing state, cleanup | integration tests with temp dirs | Gap | Must validate deterministic reload |

## Acceptance Criteria
- One Julia source file exists for every mapped Python source file, and one Julia example exists for every Python example.
- All examples execute in Julia without runtime errors on CPU.
- Parity thresholds against Python baseline:
  - use combined relative + absolute metrics
  - use denominator-floored relative metrics in deep-null/high-contrast regimes
  - use case-specific documented overrides where needed
  - threshold source-of-truth is machine-readable config in `test/parity/thresholds/` and documented in `docs/parity_thresholds.md`
- Core kernels (`prop_propagate`, `prop_lens`, `prop_dm`, interpolation path) pass allocation and type-stability checks in benchmark tests.
- Hot-path kernels meet dispatch/inference gate: no unintended dynamic dispatch or `Any` inference in designated inner-loop paths.
- Comprehensive benchmark suite is present and produces reproducible Python-vs-Julia steady-state reports.
- Julia TTFx is measured separately from steady-state runtime and reported independently.

## Benchmarking Plan (Python vs Julia)
Decision reference: `D-0029` in `docs/compat_decisions.md`.

### 1. Benchmark Categories
- Microbenchmarks (kernel-level):
  - `prop_shift_center`, `prop_qphase`, `prop_ptp`, `prop_wts`, `prop_stw`, `prop_lens`
  - interpolation/mapping kernels (`prop_magnify`, `prop_szoom`, `prop_resamplemap`, `prop_rotate`)
  - map generation (`prop_psd_errormap`) with seeded RNG
- Mesobenchmarks (workflow-level):
  - representative wavefront pipelines combining multiple kernels
  - DM and map workflows
- End-to-end benchmarks (user-level):
  - selected examples (`simple_prescription`, `simple_telescope`, `psdtest`, `run_coronagraph_dm`, `multi_example`)

### 2. Workload Matrix
- Grid sizes: `256`, `512`, `1024` (extend as needed)
- Baseline: patched Python 3.3.4 executable reference
- Backends:
  - CPU required
  - GPU optional/nightly (where implementation exists)
- Thread counts:
  - controlled, fixed for both Python and Julia runs

### 3. Measurement Protocol
- Hardware/environment controls:
  - pin thread counts and affinity where possible
  - record CPU/GPU model, Julia/Python versions, dependency versions
  - disable opportunistic thermal/power scaling where feasible on dedicated benchmark hosts
- Python timing:
  - warm up once before timed repeats
  - use repeated timing with robust summary statistics
- Julia timing (steady-state):
  - execute warmup calls to force compilation before measurement
  - run measured trials in persistent Julia process
  - use benchmark tooling for robust statistics
- Julia TTFx timing (separate track):
  - measure cold-start and first-call compile latency independently
  - do not combine with steady-state runtime comparisons
- Statistical policy:
  - use fixed sample-count minimums per case (default `>= 30` timed samples)
  - report median and p90 as primary statistics; treat minimum as informational only
  - keep outlier handling explicit and consistent across Python and Julia

### 4. TTFx Exclusion Rule (Required)
- Primary Python-vs-Julia performance claims are based on steady-state timings only.
- Any chart/table comparing Python and Julia must explicitly state that Julia TTFx is excluded.
- TTFx results are published in a separate section labeled `Cold Start / TTFx`.

### 5. Precompile Optimization Workflow
- Use representative workloads to generate/maintain precompile statements.
- Prefer Julia-native startup behavior by default; add explicit precompile directives only if TTFx data justifies them.
- Re-run TTFx suite after precompile updates.
- Accept precompile changes when:
  - TTFx improves materially, and
  - steady-state runtime is not regressed.
- Prefer codified precompile workloads over ad hoc statements so updates remain reviewable and reproducible.

### 6. Reported Metrics
- Steady-state:
  - median, p90, and best runtime
  - speedup ratio (Python / Julia)
  - allocations and allocated bytes (Julia)
- TTFx:
  - startup latency
  - first-call latency by benchmark case

### 7. Artifacts And Reproducibility
- Store benchmark config + metadata with each run:
  - commit hash, decision set, baseline tag, backend, hardware, thread count
- Emit machine-readable outputs (`JSON`/`CSV`) and human summary (`Markdown`).
- Keep large/raw artifacts outside default git history unless explicitly approved.
- Include explicit run tags:
  - `steady_state` for Python-vs-Julia runtime comparison
  - `cold_start` for Julia startup/TTFx diagnostics

### 8. Benchmark Gates
- No regression gate (steady-state): fail if core benchmark median regresses beyond agreed threshold (default `10%`) unless waived.
- Comparison target:
  - maintain or improve versus Python baseline per `D-0013` target.
- TTFx gate:
  - track trend over time; precompile work should reduce or hold TTFx, not increase unexpectedly.

### 9. Benchmark Harness Layout
- `bench/python/`:
  - baseline runners for micro, meso, and e2e workloads
  - emits JSON/CSV artifacts with workload metadata
- `bench/julia/steady_state/`:
  - persistent-process benchmark runners for fair runtime comparison
  - includes warmup path before all timed measurements
- `bench/julia/cold_start/`:
  - fresh-process scripts for startup/TTFx measurement only
  - never merged into Python-vs-Julia runtime charts
- `bench/common/`:
  - workload definitions, seed management, and shared metadata schema
- `bench/reports/`:
  - generated markdown summaries and comparison tables

### 10. Execution Policy
- Python-vs-Julia comparison command must run steady-state suites only.
- Cold-start/TTFx command must run Julia cold-start suite only.
- CI benchmark jobs must publish both artifacts, but the performance pass/fail gate uses steady-state artifacts.
- TTFx artifacts are required for trend tracking and precompile validation, not for Python-runtime parity claims.

## Suggested Julia Dependencies
- Core: `LinearAlgebra`, `Statistics`, `Random`, `FFTW`, `AbstractFFTs`
- IO (required): `FITSIO` (`FITSIO.jl`)
- Plotting (default for examples): `Plots` (`Plots.jl`)
- Optional backend/perf: `CUDA`, `KernelAbstractions`, `AcceleratedKernels`, `MKL`
- Testing/benchmarking: `Test`, `BenchmarkTools`

## One-to-One Filename Map
Note: mapping below guarantees file-level traceability; implementation may and should diverge internally to use Julia dispatch/traits.
### Core package files (`proper/*`)
| Python file | Julia file | Port notes |
| --- | --- | --- |
| `proper/.use_ffti` | `data/.use_ffti` | Compatibility marker; reassess necessity |
| `proper/__init__.py` | `src/Proper.jl` | Package entrypoint and exports |
| `proper/cubic_conv_c.c` | — | Historical upstream C reference; Julia runtime uses `src/libcconv.jl` |
| `proper/cubic_conv_threaded_c.c` | — | Historical upstream C reference; Julia runtime uses `src/libcconvthread.jl` |
| `proper/influence_dm5v2_1.fits` | `data/influence_dm5v2_1.fits` | Reference data asset |
| `proper/libcconv.py` | `src/libcconv.jl` | Direct function/module port |
| `proper/libcconvthread.py` | `src/libcconvthread.jl` | Direct function/module port |
| `proper/libszoom.py` | `src/libszoom.jl` | Direct function/module port |
| `proper/prop_8th_order_mask.py` | `src/prop_8th_order_mask.jl` | Direct function/module port |
| `proper/prop_add_phase.py` | `src/prop_add_phase.jl` | Direct function/module port |
| `proper/prop_add_wavefront.py` | `src/prop_add_wavefront.jl` | Direct function/module port |
| `proper/prop_begin.py` | `src/prop_begin.jl` | Direct function/module port |
| `proper/prop_circular_aperture.py` | `src/prop_circular_aperture.jl` | Direct function/module port |
| `proper/prop_circular_obscuration.py` | `src/prop_circular_obscuration.jl` | Direct function/module port |
| `proper/prop_compile_c.py` | `src/prop_compile_c.jl` | Direct function/module port |
| `proper/prop_cubic_conv.py` | `src/prop_cubic_conv.jl` | Direct function/module port |
| `proper/prop_define_entrance.py` | `src/prop_define_entrance.jl` | Direct function/module port |
| `proper/prop_dftidefs.py` | `src/prop_dftidefs.jl` | Direct function/module port |
| `proper/prop_divide.py` | `src/prop_divide.jl` | Direct function/module port |
| `proper/prop_dm.py` | `src/prop_dm.jl` | Direct function/module port |
| `proper/prop_ellipse.py` | `src/prop_ellipse.jl` | Direct function/module port |
| `proper/prop_elliptical_aperture.py` | `src/prop_elliptical_aperture.jl` | Direct function/module port |
| `proper/prop_elliptical_obscuration.py` | `src/prop_elliptical_obscuration.jl` | Direct function/module port |
| `proper/prop_end.py` | `src/prop_end.jl` | Direct function/module port |
| `proper/prop_end_savestate.py` | `src/prop_end_savestate.jl` | Direct function/module port |
| `proper/prop_errormap.py` | `src/prop_errormap.jl` | Direct function/module port |
| `proper/prop_execute_multi.py` | `src/prop_execute_multi.jl` | Direct function/module port |
| `proper/prop_ffti.py` | `src/prop_ffti.jl` | Direct function/module port |
| `proper/prop_fftw.py` | `src/prop_fftw.jl` | Direct function/module port |
| `proper/prop_fftw_wisdom.py` | `src/prop_fftw_wisdom.jl` | Direct function/module port |
| `proper/prop_fit_dm.py` | `src/prop_fit_dm.jl` | Direct function/module port |
| `proper/prop_fit_zernikes.py` | `src/prop_fit_zernikes.jl` | Direct function/module port |
| `proper/prop_fits_read.py` | `src/prop_fits_read.jl` | Direct function/module port |
| `proper/prop_fits_write.py` | `src/prop_fits_write.jl` | Direct function/module port |
| `proper/prop_get_amplitude.py` | `src/prop_get_amplitude.jl` | Direct function/module port |
| `proper/prop_get_beamradius.py` | `src/prop_get_beamradius.jl` | Direct function/module port |
| `proper/prop_get_distancetofocus.py` | `src/prop_get_distancetofocus.jl` | Direct function/module port |
| `proper/prop_get_fratio.py` | `src/prop_get_fratio.jl` | Direct function/module port |
| `proper/prop_get_gridsize.py` | `src/prop_get_gridsize.jl` | Direct function/module port |
| `proper/prop_get_nyquistsampling.py` | `src/prop_get_nyquistsampling.jl` | Direct function/module port |
| `proper/prop_get_phase.py` | `src/prop_get_phase.jl` | Direct function/module port |
| `proper/prop_get_refradius.py` | `src/prop_get_refradius.jl` | Direct function/module port |
| `proper/prop_get_sampling.py` | `src/prop_get_sampling.jl` | Direct function/module port |
| `proper/prop_get_sampling_arcsec.py` | `src/prop_get_sampling_arcsec.jl` | Direct function/module port |
| `proper/prop_get_sampling_radians.py` | `src/prop_get_sampling_radians.jl` | Direct function/module port |
| `proper/prop_get_wavefront.py` | `src/prop_get_wavefront.jl` | Direct function/module port |
| `proper/prop_get_wavelength.py` | `src/prop_get_wavelength.jl` | Direct function/module port |
| `proper/prop_get_z.py` | `src/prop_get_z.jl` | Direct function/module port |
| `proper/prop_hex_wavefront.py` | `src/prop_hex_wavefront.jl` | Direct function/module port |
| `proper/prop_hex_zernikes.py` | `src/prop_hex_zernikes.jl` | Direct function/module port |
| `proper/prop_init_savestate.py` | `src/prop_init_savestate.jl` | Direct function/module port |
| `proper/prop_irregular_polygon.py` | `src/prop_irregular_polygon.jl` | Direct function/module port |
| `proper/prop_is_statesaved.py` | `src/prop_is_statesaved.jl` | Direct function/module port |
| `proper/prop_lens.py` | `src/prop_lens.jl` | Direct function/module port |
| `proper/prop_load_fftw_wisdom.py` | `src/prop_load_fftw_wisdom.jl` | Direct function/module port |
| `proper/prop_magnify.py` | `src/prop_magnify.jl` | Direct function/module port |
| `proper/prop_multiply.py` | `src/prop_multiply.jl` | Direct function/module port |
| `proper/prop_noll_zernikes.py` | `src/prop_noll_zernikes.jl` | Direct function/module port |
| `proper/prop_pixellate.py` | `src/prop_pixellate.jl` | Direct function/module port |
| `proper/prop_polygon.py` | `src/prop_polygon.jl` | Direct function/module port |
| `proper/prop_print_zernikes.py` | `src/prop_print_zernikes.jl` | Direct function/module port |
| `proper/prop_propagate.py` | `src/prop_propagate.jl` | Direct function/module port |
| `proper/prop_psd_errormap.py` | `src/prop_psd_errormap.jl` | Direct function/module port |
| `proper/prop_ptp.py` | `src/prop_ptp.jl` | Direct function/module port |
| `proper/prop_qphase.py` | `src/prop_qphase.jl` | Direct function/module port |
| `proper/prop_radius.py` | `src/prop_radius.jl` | Direct function/module port |
| `proper/prop_readmap.py` | `src/prop_readmap.jl` | Direct function/module port |
| `proper/prop_rectangle.py` | `src/prop_rectangle.jl` | Direct function/module port |
| `proper/prop_rectangular_aperture.py` | `src/prop_rectangular_aperture.jl` | Direct function/module port |
| `proper/prop_rectangular_obscuration.py` | `src/prop_rectangular_obscuration.jl` | Direct function/module port |
| `proper/prop_resamplemap.py` | `src/prop_resamplemap.jl` | Direct function/module port |
| `proper/prop_rotate.py` | `src/prop_rotate.jl` | Direct function/module port |
| `proper/prop_rounded_rectangle.py` | `src/prop_rounded_rectangle.jl` | Direct function/module port |
| `proper/prop_run.py` | `src/prop_run.jl` | Direct function/module port |
| `proper/prop_run_multi.py` | `src/prop_run_multi.jl` | Direct function/module port |
| `proper/prop_savestate.py` | `src/prop_savestate.jl` | Direct function/module port |
| `proper/prop_select_propagator.py` | `src/prop_select_propagator.jl` | Direct function/module port |
| `proper/prop_set_antialiasing.py` | `src/prop_set_antialiasing.jl` | Direct function/module port |
| `proper/prop_shift_center.py` | `src/prop_shift_center.jl` | Direct function/module port |
| `proper/prop_sinc.py` | `src/prop_sinc.jl` | Direct function/module port |
| `proper/prop_state.py` | `src/prop_state.jl` | Direct function/module port |
| `proper/prop_stw.py` | `src/prop_stw.jl` | Direct function/module port |
| `proper/prop_szoom.py` | `src/prop_szoom.jl` | Direct function/module port |
| `proper/prop_szoom_c.c` | — | Historical upstream C reference; Julia runtime uses `src/prop_szoom.jl` |
| `proper/prop_table.py` | `src/prop_table.jl` | Direct function/module port |
| `proper/prop_use_ffti.py` | `src/prop_use_ffti.jl` | Direct function/module port |
| `proper/prop_use_fftw.py` | `src/prop_use_fftw.jl` | Direct function/module port |
| `proper/prop_wavefront.py` | `src/prop_wavefront.jl` | Direct function/module port |
| `proper/prop_writemap.py` | `src/prop_writemap.jl` | Direct function/module port |
| `proper/prop_wts.py` | `src/prop_wts.jl` | Direct function/module port |
| `proper/prop_zernikes.py` | `src/prop_zernikes.jl` | Direct function/module port |
| `proper/switch_set.py` | `src/switch_set.jl` | Direct function/module port |

### Example files (`proper/examples/*`)
| Python file | Julia file | Port notes |
| --- | --- | --- |
| `proper/examples/__init__.py` | `examples/__init__.jl` | Example loader/namespace parity |
| `proper/examples/coronagraph.py` | `examples/coronagraph.jl` | Direct example prescription/demo port |
| `proper/examples/coronagraph_demo.py` | `examples/coronagraph_demo.jl` | Direct example prescription/demo port |
| `proper/examples/example_system.py` | `examples/example_system.jl` | Direct example prescription/demo port |
| `proper/examples/hubble_simple.py` | `examples/hubble_simple.jl` | Direct example prescription/demo port |
| `proper/examples/microscope.py` | `examples/microscope.jl` | Direct example prescription/demo port |
| `proper/examples/multi_example.py` | `examples/multi_example.jl` | Direct example prescription/demo port |
| `proper/examples/occulter_demo.py` | `examples/occulter_demo.jl` | Direct example prescription/demo port |
| `proper/examples/psdtest.py` | `examples/psdtest.jl` | Direct example prescription/demo port |
| `proper/examples/run_coronagraph.py` | `examples/run_coronagraph.jl` | Direct example prescription/demo port |
| `proper/examples/run_coronagraph_dm.py` | `examples/run_coronagraph_dm.jl` | Direct example prescription/demo port |
| `proper/examples/run_example.py` | `examples/run_example.jl` | Direct example prescription/demo port |
| `proper/examples/run_occulter.py` | `examples/run_occulter.jl` | Direct example prescription/demo port |
| `proper/examples/simple_prescription.py` | `examples/simple_prescription.jl` | Direct example prescription/demo port |
| `proper/examples/simple_telescope.py` | `examples/simple_telescope.jl` | Direct example prescription/demo port |
| `proper/examples/talbot.py` | `examples/talbot.jl` | Direct example prescription/demo port |
| `proper/examples/talbot_correct.py` | `examples/talbot_correct.jl` | Direct example prescription/demo port |
| `proper/examples/talbot_correct_demo.py` | `examples/talbot_correct_demo.jl` | Direct example prescription/demo port |
| `proper/examples/talbot_demo.py` | `examples/talbot_demo.jl` | Direct example prescription/demo port |
| `proper/examples/telescope.py` | `examples/telescope.jl` | Direct example prescription/demo port |
| `proper/examples/telescope_dm.py` | `examples/telescope_dm.jl` | Direct example prescription/demo port |
| `proper/examples/testmulti1.py` | `examples/testmulti1.jl` | Direct example prescription/demo port |
| `proper/examples/testmulti2.py` | `examples/testmulti2.jl` | Direct example prescription/demo port |

## Example Port Checklist
- [ ] `coronagraph.jl`
- [ ] `coronagraph_demo.jl`
- [ ] `example_system.jl`
- [ ] `hubble_simple.jl`
- [ ] `microscope.jl`
- [ ] `multi_example.jl`
- [ ] `occulter_demo.jl`
- [ ] `psdtest.jl`
- [ ] `run_coronagraph.jl`
- [ ] `run_coronagraph_dm.jl`
- [ ] `run_example.jl`
- [ ] `run_occulter.jl`
- [ ] `simple_prescription.jl`
- [ ] `simple_telescope.jl`
- [ ] `talbot.jl`
- [ ] `talbot_correct.jl`
- [ ] `talbot_correct_demo.jl`
- [ ] `talbot_demo.jl`
- [ ] `telescope.jl`
- [ ] `telescope_dm.jl`
- [ ] `testmulti1.jl`
- [ ] `testmulti2.jl`

## Execution Order Recommendation
1. Complete Phase 0 preflight decisions and freeze contracts (API, numerics, backend traits, parity harness, CI).
2. Port foundational wavefront + propagation core and run `simple_prescription` parity.
3. Port aperture/mask/geometry and run telescope + microscope parity.
4. Port map/PSD/interpolation stack and run `psdtest` parity.
5. Port coronagraph and DM stack and run coronagraph + DM + multi-run examples.
6. Enable and validate GPU backend for selected core examples, then expand to full suite.
7. Execute full parity closure pass (replace all fallback physics paths and close Python parity gaps).
8. Perform MATLAB/manual semantic reconciliation and finalize release validation.
