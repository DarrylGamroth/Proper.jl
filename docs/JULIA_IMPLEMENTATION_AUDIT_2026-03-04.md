# Proper.jl Full Implementation Audit (2026-03-04)

## Scope
Audit goals requested:
- Identify opportunities to improve type selectors, multiple dispatch, and traits.
- Identify idiomatic Julia opportunities (array indexing, comprehensions/`map`, tuples, etc.).
- Find hard-coded concrete numeric types (especially `Float64`) and non-parametric patterns.
- Identify allocation hotspots and ways to reduce allocations.
- Review potential use of `Bumper.jl`, `AllocArrays.jl`, `KernelAbstractions.jl`, and `AcceleratedKernels.jl`.

Audited paths:
- `src/**/*.jl`
- relevant design docs (`docs/backend_traits.md`, `docs/compat_decisions.md`)

## Method
1. Static scan for dynamic-dispatch and non-parametric patterns (`kwargs...`, `Any`, hard-coded `Float64`, `collect`, `Dict{...,Any}`, `Vector{Any}`).
2. Manual code review of core and hotspot modules.
3. Targeted inference and allocation probes (`@inferred`, `@code_warntype`, `@allocated`) on representative kernels.

## Executive Summary
The codebase has solid early structure (typed `WaveFront`, compatibility policy types, trait scaffolding), but runtime behavior is still largely dynamic in the performance-critical paths.

Top issues:
- Trait/policy architecture is defined but minimally used by runtime kernels.
- Several functions still hard-code `Float64`/CPU `Array` construction, breaking backend preservation and parametric design.
- `prop_fits_read`/`prop_readmap` are type-unstable (`Any` inference).
- High allocation pressure remains in geometry/mask generation and PSD/hex workflows.

This is recoverable with a focused refactor that keeps public `prop_*` APIs stable while moving internals to typed options + mutating kernels.

## High-Priority Findings

### F1. Trait/policy architecture is not yet driving runtime kernels
Severity: High

Evidence:
- `RunContext` exists but is not threaded through `prop_*` runtime paths: `src/core/context.jl:3-30`.
- Traits are placeholders with generic defaults only: `src/core/traits.jl:1-23`.
- Most kernel entry points dispatch only on `WaveFront` + dynamic kwargs, not policy/backend traits.

Impact:
- Accepted decisions (`D-0017`, `D-0019`, `D-0028`) are only partially realized.
- GPU/backend specialization cannot be introduced cleanly without reworking call chains.

Recommendation:
- Add internal mutating kernels with typed context:
  - `prop_resamplemap!(out, wf, dmap, opts, ctx)`
  - `prop_rotate!(out, image, θ, opts, ctx)`
  - `prop_magnify!(out, image, mag, opts, ctx)`
- Keep `prop_*` wrappers for compatibility, but normalize kwargs once and call typed kernels.

---

### F2. Type-instability in FITS/readmap path
Severity: High

Evidence:
- `prop_fits_read` initializes `arr = nothing` then writes in closure and returns `Float64.(arr)`:
  - `src/prop_fits_read.jl:4-18`
- `@code_warntype` shows `Body::ANY` for `prop_fits_read` and `prop_readmap`.

Impact:
- Inference barrier on map I/O path.
- Avoidable dynamic dispatch and potential downstream instability.

Recommendation:
- Return concretely typed results directly from the FITS `do` block.
- Split methods:
  - `prop_fits_read(fname)` -> `AbstractArray{T}`
  - `prop_fits_read_with_header(fname)` -> tuple with typed header wrapper
- Avoid forcing `Float64`; preserve element type or accept explicit `T` conversion argument.

---

### F3. Hard-coded `Float64` and CPU arrays break parametric/GPU intent
Severity: High

Evidence:
- `src/core/context.jl:12-17` defaults to `Array{Float64,2}`.
- `src/prop_fit_zernikes.jl:3,29-31,50,74` uses `Matrix{Float64}`, `Float64[]`, `zeros(Float64,...)`, `ones(Float64,...)`.
- `src/libszoom.jl:12-13` uses `Vector{Float64}`.
- `src/prop_rounded_rectangle.jl:12` uses `zeros(Float64,...)`.
- `src/prop_fits_read.jl:17` hard-casts FITS data to `Float64`.

Impact:
- Violates backend-preserving policy and precision genericity.
- Blocks clean GPU dispatch where values should track input eltype/backend.

Recommendation:
- Replace with `similar`, `similar_type`, or `promote_type`-driven allocations.
- Use `RT = real(eltype(wf.field))` (or propagated `T`) consistently.
- Reserve explicit `Float64` only for deliberately fixed-format artifacts.

---

### F4. Dynamic containers (`Any`) in multirun path
Severity: High

Evidence:
- `src/prop_run_multi.jl:7` uses `Vector{Any}` for outputs.
- Array stack eltype discovered at runtime from `outputs[1]`: `src/prop_run_multi.jl:18-21`.

Impact:
- Type instability and dynamic dispatch in a core orchestration path.
- Harder to optimize threaded execution.

Recommendation:
- First-pass evaluate one run to determine output type and dimensions, then allocate typed output container.
- Or require homogeneous output type contract and enforce with typed buffering.

## Medium-Priority Findings

### F5. Symbol/string propagation state machine in hot propagation path
Severity: Medium

Evidence:
- `WaveFront` stores `reference_surface`, `beam_type_old`, `propagator_type` as `Symbol`: `src/core/types.jl:22-24`.
- Runtime string-to-symbol transitions: `src/prop_lens.jl:40`, `src/prop_select_propagator.jl:6`.

Impact:
- Avoidable runtime symbol/string work.
- Prevents richer compile-time dispatch on state transitions.

Recommendation:
- Replace `Symbol` state fields with compact `@enum` fields.
- Dispatch transition math by state enum pair.
- Keep user-facing names unchanged.

---

### F6. Extensive kwargs/haskey dynamic option parsing
Severity: Medium

Evidence:
- 68 `kwargs...` definitions in `src/`.
- Repeated runtime `haskey`/symbol normalization patterns in many kernels:
  - `src/switch_set.jl:2-10`
  - `src/prop_psd_errormap.jl:3-107`
  - `src/prop_hex_wavefront.jl:16-22`
  - `src/prop_rectangle.jl:8-10`

Impact:
- Dynamic branching overhead and duplicated option parsing logic.
- Hard to ensure type-stable hot kernels.

Recommendation:
- One-time normalize kwargs at boundary into typed options structs (e.g., `ResampleOptions`, `MaskOptions`, `PSDOptions`).
- Dispatch kernels on option type + policy trait.

---

### F7. Return-type unions from mode/option booleans
Severity: Medium

Evidence:
- `prop_fit_zernikes(...; fit=false)` returns either `Vector` or `(Vector, Matrix)`:
  - `src/prop_fit_zernikes.jl:49-69`
- `prop_cubic_conv(...; grid::Bool=...)` uses runtime bool to choose topology:
  - `src/prop_cubic_conv.jl:54-60`

Impact:
- Inference friction and less predictable call sites.

Recommendation:
- Split into separate method names or dispatch on topology/type selectors instead of runtime bool flags.

---

### F8. Avoidable temporary allocations and non-idiomatic axis construction
Severity: Medium

Evidence:
- Frequent `collect(...)` in axis construction and option conversion (`13` sites), e.g.:
  - `src/core/grid.jl:1-4,19-22`
  - `src/prop_polygon.jl:14`
  - `src/prop_psd_errormap.jl:27`
- Broadcast temporaries in geometry kernels (e.g., ellipse):
  - `src/prop_ellipse.jl:56-59`

Impact:
- Extra allocation pressure and cache traffic.

Recommendation:
- Prefer explicit loops for hotspot math.
- Use typed preallocated vectors for axes.
- Replace `collect`-heavy expression chains with in-place fill kernels.

## Low-Priority Idiomatic Opportunities

### F9. Memory access order refinements
Severity: Low

Evidence:
- Some loops write arrays with non-optimal index order in column-major storage:
  - `src/prop_ellipse.jl:48-53` writes `x_local[j, i]` / `y_local[j, i]` with `i` as inner index.

Recommendation:
- Standardize loop order to iterate first index (`i`) in inner loop for contiguous writes where possible.

---

### F10. Additional dispatch opportunities
Severity: Low

Evidence:
- `prop_zernikes` and `prop_hex_zernikes` normalize scalar/vector inputs at runtime with `isa`/`collect`:
  - `src/prop_zernikes.jl:149-150`
  - `src/prop_hex_zernikes.jl:76-77`

Recommendation:
- Add separate scalar/vector methods to avoid runtime normalization branches.

## Allocation Probe Results (Representative)

Measured using warmed-up `@allocated` calls on representative inputs.

| Function | Approx. allocated bytes |
| --- | ---: |
| `prop_radius` | 533,064 |
| `prop_magnify` | 754,136 |
| `prop_rotate` | 524,376 |
| `prop_resamplemap` | 524,376 |
| `prop_cubic_conv(a, y, x)` | 0 |
| `prop_cubic_conv(...; grid=true)` | 2,361,672 |
| `prop_ellipse` | 562,928 |
| `prop_rectangle` | 524,920 |
| `prop_irregular_polygon` | 16,600,464 |
| `prop_zernikes(wf, 8)` | 4,203,064 |
| `prop_zernikes` apply variant | 102,295,520 |
| `prop_fit_zernikes` | 50,356,384 |
| `prop_psd_errormap(...; NO_APPLY=true)` | 320,736,776 |
| `prop_8th_order_mask(...; circular=true)` | 124,284,608 |
| `prop_hex_wavefront(...; no_apply=true)` | 73,194,840 |

Interpretation:
- Hot kernels are not close to zero-allocation steady-state targets from `D-0013`/`D-0021`.
- Most large allocations are from temporary array materialization and repeated full-frame intermediates.

## Package Review: Bumper.jl and AllocArrays.jl

### Bumper.jl
Fit:
- Good for short-lived scratch arrays in tight CPU loops via task-local slab allocation.
- Best when allocations are numerous and bounded within a non-escaping region.

Where it could help here:
- Temporary buffers in `prop_ellipse`, `prop_psd_errormap`, `prop_hex_wavefront`, and Zernike assembly paths.

Limits:
- Not a substitute for structural refactors (`kwargs` normalization, dispatch cleanup).
- Not a GPU strategy.
- Unsuitable for buffers that must escape call scope.

Recommendation:
- Use only after kernel interfaces are converted to mutating typed forms.
- Keep optional behind internal allocator trait; do not make behavior depend on Bumper.

### AllocArrays.jl
Fit:
- Arena-style allocation can reduce GC pressure for repeated temporary arrays.

Where it could help:
- Repeated temporary matrices in PSD/map synthesis and segmentation loops.

Limits:
- Similar to Bumper: CPU memory-management optimization, not a dispatch/backend architecture solution.
- Additional dependency complexity should be justified by benchmark deltas.

Recommendation:
- Run a scoped benchmark spike (PSD + hex workloads) before adopting.
- Keep optional and internal; avoid exposing allocator semantics in public API.

## Package Review: KernelAbstractions.jl / AcceleratedKernels.jl

### KernelAbstractions.jl (KA)
Strong candidate for explicit backend-portable kernels in:
- `prop_resamplemap`, `prop_magnify`, `prop_rotate`
- `prop_cubic_conv` grid sampling
- geometric mask fills (`prop_rectangle`, `prop_ellipse`, irregular polygon sampling)

Approach:
- Introduce `*_kernel!` methods and dispatch via backend traits.
- Keep CPU fallback loops for parity/debug.

### AcceleratedKernels.jl (AK)
Good candidate for high-level fused map/reduction operations where API coverage matches needs:
- reduction/scaling sub-steps in PSD normalization and map postprocessing
- replacing broadcast chains that currently materialize intermediates

Approach:
- Use AK selectively where it simplifies code and benchmarks positively.
- Do not force AK for complex gather-heavy interpolation where KA custom kernels are clearer.

## Recommended Refactor Plan (Audit-Driven)

### Stage A: Type and inference hardening
1. Fix `prop_fits_read`/`prop_readmap` inference (`Any` elimination).
2. Remove hard-coded `Float64` and `ComplexF64` from core math paths.
3. Split union-return methods (`prop_fit_zernikes`, boolean mode switches) into dispatch-specialized methods.

### Stage B: Option typing and dispatch cleanup
1. Introduce typed options structs per kernel family.
2. Normalize kwargs once in wrappers; move kernels to typed signatures.
3. Replace symbol/string propagation state fields with enums.

### Stage C: Allocation control
1. Convert hotspot kernels to mutating `*_!` APIs with reusable workspace.
2. Add allocation regression tests for selected kernels.
3. Evaluate optional Bumper/AllocArrays acceleration after structural cleanup.

### Stage D: Backend acceleration
1. Implement trait-routed CPU + KA kernels for resample/rotate/magnify/cubic conv.
2. Add AK where it removes large intermediate arrays cleanly.
3. Add GPU tests with `CUDA.allowscalar(false)` for implemented paths.

## Immediate Actionable Checklist
- [ ] Eliminate `Any` in FITS/readmap path.
- [ ] Remove all hard-coded `Float64` array allocations in hot kernels.
- [ ] Replace `Vector{Any}` in `prop_run_multi`.
- [ ] Introduce first typed options struct and one `*_!` kernel path (`prop_resamplemap!`).
- [ ] Add CI checks for inference (`@inferred`) and allocation budget on selected kernels.

