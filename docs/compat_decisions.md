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
- Consequences:
  - Automated golden-data generation and CI parity checks run against Python outputs.
  - Divergences from Python must be documented in this file with justification and tests.
  - Behavior changes are tracked by decisions/tests rather than runtime compatibility mode flags.

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
- Status: Superseded by D-0035
- Decision:
  - Historical: default mode `:python334`, optional `:corrected`.
  - Current: runtime compatibility modes removed; use patched baseline parity policy (D-0035).

## D-0007: Numerical Convention Contract
- Date: 2026-03-04
- Status: Accepted
- Decision:
  - Match Python 3.3.4 normalization and centering behavior.
  - Freeze coordinate/pixel-center conventions per PROPER behavior and test them explicitly (even/odd sizes).
  - Expose all numerically meaningful deviations via explicit decisions and tests.

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
- Status: Superseded by D-0035
- Context:
  - Repeated `compat_mode` branching inside inner loops adds overhead and increases behavioral complexity.
  - Mode decisions should be deterministic and easy to audit.
- Decision:
  - Historical: accepted mode exactly once in constructor and passed typed policy internally.
  - Current: compatibility modes removed; no mode dispatch remains.
- Consequences:
  - Public `prop_*` APIs and context constructors have no compatibility-mode parameters.
  - Behavior differences are represented by explicit code paths and decision-log entries.

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
  - Introduce a typed run context carrying backend traits, RNG, and workspace handles.
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
- Status: Superseded by D-0035
- Context:
  - Historical policy before runtime compatibility mode removal.
- Decision:
  - Historical: `:corrected` tracked as behavior-equivalent to `:python334`.
- Consequences:
  - Replaced by single baseline policy in D-0035.

## D-0033: Phase 9 Hotspot Reconciliation Outcomes
- Date: 2026-03-04
- Status: Accepted
- Context:
  - Phase 9 requires explicit reconciliation of known Python translation defects using MATLAB/manual semantics where justified.
- Decision:
  - `prop_resamplemap`: keep independent `xshift`/`yshift` semantics.
  - `prop_end`: keep integer-safe centered extraction semantics (do not preserve Python float-slice defect).
  - `prop_state`: restore full wavefront state fields in-place (do not preserve Python local-rebind defect).
  - `prop_psd_errormap`: keep explicit Julia FFT path without backend-toggle side effects; retain documented parity quirks only where needed for executable parity.
- Consequences:
  - Reconciled behaviors are test-covered in `test/test_phase9_semantic_reconciliation.jl`.
  - Future changes to these hotspots require decision-log updates and parity evidence.

## D-0034: Parity Baseline Patch For `prop_resamplemap`
- Date: 2026-03-04
- Status: Accepted
- Context:
  - Python 3.3.4 source included a `prop_resamplemap` Y-shift defect (`y += yc - xshift/pixscale`) that disagreed with MATLAB and intended semantics.
  - Project direction is to match upstream behavior after applying this known fix, not preserve the historical defect.
- Decision:
  - Patch the local Python parity baseline to use `y += yc - yshift/pixscale`.
  - Regenerate parity artifacts from the patched baseline and treat that output as canonical executable reference in this repository.
  - Keep Julia `prop_resamplemap` on independent `xshift`/`yshift` semantics.
- Consequences:
  - Removes ambiguity around preserving the historical Y-shift coupling defect.
  - Parity comparisons stay executable and reproducible while aligned with intended semantics.

## D-0035: Remove Runtime Compatibility Modes
- Date: 2026-03-04
- Status: Accepted
- Context:
  - `compat_mode` no longer provided real behavioral coverage and added maintenance/documentation overhead.
  - Project direction is to keep one executable baseline and track intentional divergences explicitly.
- Decision:
  - Remove runtime `compat_mode` support from context/config/public APIs.
  - Remove compatibility policy types and mode-based dispatch plumbing.
  - Keep parity anchored to the patched Python 3.3.4 executable baseline.
- Consequences:
  - Supersedes D-0006, D-0017, and D-0032.
  - Future behavior changes require decision-log entries + tests instead of mode branches.

## D-0036: `prop_dm`/`prop_cubic_conv` Indexing Reconciliation (MATLAB/Python)
- Date: 2026-03-04
- Status: Accepted
- Context:
  - `prop_dm` was materially diverging in DM-corrected coronagraph behavior despite passing coarse parity summaries.
  - MATLAB reference (`prop_dm.m`) uses 1-based interpolation coordinates (`xc_inf = floor(nx/2)+1`, actuator grid origins with `+1` index offsets).
  - Python translation maps those coordinates into 0-based cubic-convolution coordinates.
- Decision:
  - Align `libcconv` with upstream `cubic_conv_c.c` semantics (rounding, support window, edge clamping).
  - In Julia `prop_dm`, keep interpolation-coordinate math in upstream 0-based space while converting only array indexing operations to 1-based.
  - Apply DM phase using transposed internal map layout to match current wavefront-array orientation in this codebase.
- Consequences:
  - `run_coronagraph_dm` with deterministic map input now matches Python 3.3.4 executable outputs for no-error, error-only, and DM-corrected cases.
  - Added regression coverage to prevent reintroduction of DM correction divergence.

## D-0037: WFIRST SPC Public-Data Compatibility Mapping
- Date: 2026-03-14
- Status: Accepted
- Context:
  - The executable Python WFIRST Phase B baseline in `../proper-models/wfirst_cgi/models_phaseb/python` expects legacy SPC data directories (`spc_20190130`, `spc_20181220`) that are not present in this checkout.
  - The public Roman preflight archive used by this repository contains redistributable SPC assets, but under newer 2020 directory/file names.
- Decision:
  - The local compatibility-data builder may alias public Roman SPC assets into the legacy Phase B compatibility layout used by the Python baseline.
  - Python-vs-Julia WFIRST SPC comparison harnesses in this repository compare both implementations against the same aliased compatibility root rather than the unavailable original private SPC data tree.
- Consequences:
  - SPC comparison runs remain executable and reproducible in this repository.
  - These SPC harnesses validate Julia-vs-Python parity on the shared public-data compatibility root, not historical bitwise equivalence to the unavailable original SPC datasets.

## D-0038: MATLAB As Column-Major Semantic Reference
- Date: 2026-03-16
- Status: Accepted
- Context:
  - Python 3.3.4 remains the executable parity baseline, but several recent WFIRST parity defects were caused by array-order and centering assumptions rather than optical-model logic.
  - Julia and MATLAB are both column-major, while the executable Python baseline is row-major/NumPy-oriented.
  - The repository contains both the core MATLAB PROPER reference in `../proper_v3.3.1_matlab` and the WFIRST MATLAB reference tree alongside the Python WFIRST model.
- Decision:
  - Keep Python 3.3.4 as the only executable baseline for automated parity and golden-data generation.
  - When debugging centering, transpose, cropping, indexing, or other array-order-sensitive behavior, explicitly consult the MATLAB implementations as the semantic reference for column-major correctness:
    - core PROPER: `../proper_v3.3.1_matlab`
    - WFIRST model: `../proper-models/wfirst_cgi/models_phaseb/matlab`
  - Use MATLAB selectively as a tie-breaker for semantics, not as a second runtime parity target.
- Consequences:
  - Future array-order-sensitive bugs should be evaluated against both the Python executable baseline and the MATLAB column-major reference before changing core semantics.
  - Python-order accommodations should remain local to compatibility layers or model-specific loaders unless a core PROPER routine is shown to be semantically wrong.

## D-0039: `prop_rotate` Follows Python Cubic-Interpolation Coordinates
- Date: 2026-03-16
- Status: Superseded by D-0040
- Context:
  - MATLAB `prop_rotate` delegates to `interp2` with the default one-based image-coordinate convention and centers at `fix(nx/2)+1`, `fix(ny/2)+1`.
  - The executable Python 3.3.4 baseline uses `prop_cubic_conv`, which in turn uses the upstream `cubic_conv_c.c` zero-based coordinate convention. That makes the cubic `prop_rotate` path semantically different from MATLAB on small synthetic arrays even at `theta = 0`.
  - Julia `prop_resamplemap` already matches the MATLAB source formulation, but `prop_rotate` must choose one executable convention.
- Decision:
  - Keep `prop_rotate` aligned with the executable Python baseline for both CPU and KA cubic/linear rotate kernels.
  - Document the MATLAB difference explicitly instead of trying to force both behaviors through one implementation.
  - Continue to use MATLAB as the semantic reference for column-major-sensitive routines, but not when that would break the accepted executable Python baseline without a broader compatibility decision.
- Consequences:
  - `prop_rotate` regression tests should encode the Python baseline behavior, not MATLAB `interp2` identity expectations on tiny arrays.
  - `prop_resamplemap` remains a separate MATLAB-aligned path and should not be conflated with `prop_rotate` coordinate semantics.

## D-0040: `prop_rotate` Follows MATLAB Semantics
- Date: 2026-03-16
- Status: Accepted
- Context:
  - Julia, MATLAB, and IDL are column-major, while the executable Python baseline is row-major.
  - The previous `prop_rotate` implementation mixed MATLAB-style centers with Python-style cubic-convolution sampling, which was internally inconsistent and did not cleanly match either upstream.
  - For array-order-sensitive core semantics, the project now prefers the MATLAB/IDL interpretation unless there is a stronger accepted compatibility reason not to.
- Decision:
  - `prop_rotate` follows the MATLAB `prop_rotate.m` interface and linear-path semantics:
    - default interpolation method is `linear`
    - default centers remain `fix(nx/2)+1`, `fix(ny/2)+1`
    - `MISSING` / `EXTR` controls the extrapolated fill value
    - out-of-bounds samples are not edge-clamped on the linear rotate path
  - The explicit cubic rotate path remains available and uses the upstream PROPER cubic-convolution kernel with an explicit conversion from MATLAB/Julia 1-based rotate coordinates to the kernel's 0-based coordinate space.
- Consequences:
  - `prop_rotate(a, 0)` is again an identity transform for ordinary matrices under the default settings.
  - Existing callers that relied on the previous implicit cubic default must request `METH="cubic"` or `CUBIC=true` explicitly.
  - Python executable parity is no longer the deciding baseline for the default `prop_rotate` path; MATLAB semantics are.
  - The cubic branch is now a coherent column-major translation of the upstream PROPER cubic kernel, but it is still not separately validated against executable MATLAB `interp2(..., 'cubic')`.

## D-0041: WFIRST Error-Map Public-Data Compatibility Uses Local Python-Order Loading
- Date: 2026-03-16
- Status: Accepted
- Context:
  - The executable Python WFIRST Phase B baseline applies optical error maps from legacy `wfirst_phaseb_*` FITS filenames under `data_phaseb/maps/`.
  - The public Roman preflight archive used in this repository contains redistributable error maps, but under newer `roman_phasec_*` names.
  - The first direct Julia `use_errors=1` port using that shared public-data root diverged materially while the no-error rows were already matching, which isolated the issue to error-map data interpretation rather than general branch logic.
- Decision:
  - Extend the local compatibility-data builder to alias the newer public Roman error-map FITS files into the legacy `wfirst_phaseb_*` names expected by the Python baseline.
  - Keep the array-order accommodation local to the WFIRST reference model: its error-map application path reads those FITS maps in Python order before resampling/applying them.
  - Do not change core `prop_errormap` semantics based on this WFIRST-specific parity issue alone.
- Consequences:
  - Representative `use_errors=1` HLC and SPC rows now match the executable Python baseline on the shared public-data compatibility root at the same numerical-fidelity scale as the rest of the WFIRST matrix.
  - The repository remains explicit that this is parity against the executable Python baseline on shared public data, not historical bitwise equivalence to the unavailable original private error-map tree.
  - Any future core `prop_errormap` semantic change still requires core PROPER evidence and should not be inferred from WFIRST alone.

## D-0042: WFIRST `hlc_erkin` Uses A Public-Data Compatibility Alias
- Date: 2026-03-16
- Status: Accepted
- Context:
  - The executable Python WFIRST Phase B baseline includes an `hlc_erkin` branch that expects a private legacy asset tree under `hlc_20190206_v3/` with `dsn17d_run2_pup310_fpm2048_*` filenames.
  - That original private dataset is not available in this repository, so direct historical parity against the original `hlc_erkin` files is not possible here.
  - The same repository already uses explicit public-data compatibility mappings for SPC and WFIRST error-map validation when the original private trees are unavailable.
- Decision:
  - Extend the local compatibility-data builder to synthesize the legacy `hlc_20190206_v3/` layout expected by the Python baseline from the nearest public HLC asset family available in the Roman preflight archive, `hlc_20190210b`.
  - Validate Julia-vs-Python `hlc_erkin` parity against that shared public-data compatibility root rather than the unavailable original private `hlc_20190206_v3` dataset.
  - Document this explicitly as a feature-parity and numeric-fidelity result on the shared alias, not as a claim of bitwise recovery of the private historical asset set.
- Consequences:
  - `compact_hlc_erkin` and `full_hlc_erkin` are executable and numerically matched in this repository using the same public-data compatibility root for both Python and Julia.
  - The WFIRST configuration matrix can mark `hlc_erkin` as covered for this repository's shared public-data baseline.
  - Future availability of the original private `hlc_20190206_v3` tree would justify a separate historical-parity check, but it is not required for the current public, reproducible validation surface.

## D-0043: `prop_szoom`/`prop_pixellate` Follow MATLAB Core Semantics
- Date: 2026-03-16
- Status: Accepted
- Context:
  - The current Julia `prop_szoom` implementation used a ceil/floor rounding rule inherited from the Python pure fallback, while the MATLAB reference and upstream `prop_szoom_c.c` both use nearest-integer rounding.
  - The current Julia `prop_pixellate` surface was only an integer box-downsampling helper, not the upstream PROPER PSF-integration API defined by MATLAB and Python (`prop_pixellate(image, sampling_in, sampling_out, n_out)`).
  - Julia is column-major like MATLAB/IDL, and this audit is explicitly using MATLAB as the semantic reference for core array operations unless there is a stronger accepted reason not to.
- Decision:
  - Change `prop_szoom` to use the MATLAB/C nearest-integer rounding rule.
  - Change `prop_magnify` default output sizing to use `fix` semantics for positive magnifications.
  - Add the real upstream `prop_pixellate(image, sampling_in, sampling_out, n_out=0)` overload.
  - Demote the simple integer-factor pixellation helper to an internal-only helper for tests/benchmarks instead of keeping it as a public `prop_pixellate` overload.
- Consequences:
  - Core magnification and detector-pixel integration paths now match the MATLAB formulas more closely.
  - The public `prop_pixellate` surface now matches upstream semantics instead of exposing an extra non-upstream overload.
  - The old integer-factor helper remains available only as an internal benchmark/test utility.
  - Future parity/debug work for magnification and detector-integration behavior should treat the MATLAB formulas as the semantic source unless a later accepted decision says otherwise.

## D-0044: `prop_errormap` Follows Upstream Mirror-Surface And Rotate Semantics
- Date: 2026-03-16
- Status: Accepted
- Context:
  - MATLAB and Python both apply mirror-surface maps as `exp(-4πi map / λ)`, reflecting the sign convention for reflected optical path delay.
  - MATLAB and Python both rotate error maps through the explicit cubic interpolation path after centering the map.
  - Julia had drifted on both points:
    - mirror-surface maps were applied with the wrong sign in `prop_errormap` and `prop_psd_errormap`
    - rotated error maps were using the new default linear `prop_rotate` path instead of explicit cubic semantics
- Decision:
  - Apply mirror-surface phase with negative sign in both `prop_errormap` and `prop_psd_errormap`.
  - Make the `ROTATEMAP` path in `prop_errormap` call `prop_rotate(...; METH=\"cubic\", MISSING=0.0)` explicitly.
- Consequences:
  - Core error-map application is again aligned with upstream PROPER semantics.
  - The project can keep the MATLAB-aligned default `prop_rotate` behavior while still preserving upstream `prop_errormap` behavior through an explicit cubic call.
