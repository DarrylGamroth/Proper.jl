# Backend Traits Contract

## Purpose
Define stable trait-dispatched interfaces for array backends and computational kernels.

## Decision Links
- `D-0008` GPU scope for phase 1 (accepted)
- `D-0009` backend traits and interfaces (accepted)
- `D-0012` `prop_run_multi` semantics (accepted)
- `D-0017` `compat_mode` evaluation policy (accepted)
- `D-0023` type-stable core structs (accepted)
- `D-0028` dynamic dispatch minimization (accepted)

## Status
- [x] Draft pre-filled from proposed defaults
- [x] Trait interfaces frozen
- [ ] Implementations registered

## 1. Trait Model
Proposed trait families:
- `BackendStyle(::Type{<:AbstractArray})`
- `FFTStyle(::Type{<:AbstractArray})`
- `InterpStyle(::Type{<:AbstractArray})`
- `RNGStyle(::Type)`
- `CompatPolicyStyle(::Type)` (resolved once from public `compat_mode`)

Compatibility requirements:
- Traits are dispatch aids; avoid runtime mutable global backend switches.
- Methods should preserve array backend where practical.
- Resolve compat mode at API boundary and dispatch on policy type internally.
- Hot-path methods should infer concrete types without `Any`.
- Dynamic dispatch in inner loops is disallowed unless explicitly justified.

## 2. Required Interfaces

### 2.1 FFT Interfaces
- `fft2!(A, backend_ctx)`
- `ifft2!(A, backend_ctx)`
- `fft_workspace(A, backend_ctx)` (optional)

Semantics:
- in-place where possible
- must satisfy numerics contract normalization/routine parity

### 2.2 Interpolation/Resampling Interfaces
- `resample_map!(out, in, coords, backend_ctx)`
- `magnify!(out, in, mag, backend_ctx)`
- `rotate!(out, in, angle, backend_ctx)`

Semantics:
- deterministic for fixed inputs
- behavior controlled by numerics contract, not backend-default interpolation differences

### 2.3 RNG Interfaces
- `phase_noise!(out, rng, backend_ctx)`
- Seed handling contract:
  - explicit seeded RNG for tests/parity
  - runtime default may use non-deterministic RNG

### 2.4 I/O Boundaries
- FITS I/O uses `FITSIO.jl` with explicit host/device transfers where needed.
- No direct FITS I/O from device arrays in phase 1 unless explicitly approved.

### 2.5 Compatibility Policy Interface
- `resolve_compat_policy(mode)::CompatPolicy`
- `RunContext(; compat_mode=...)` (or equivalent constructor) is the only public entry point that accepts `compat_mode`
- `CompatPolicy` propagated through workflow/context structs
- Inner loops should use policy-specialized methods, not repeated symbol comparisons

Recommended pattern:
```julia
abstract type CompatPolicy end
struct Python334Policy <: CompatPolicy end
struct CorrectedPolicy <: CompatPolicy end

resolve_compat_policy(::Val{:python334}) = Python334Policy()
resolve_compat_policy(::Val{:corrected}) = CorrectedPolicy()

# Resolve once at entry:
# ctx = RunContext(; compat_mode=:python334, ...)
# policy = ctx.policy  # already resolved in constructor
# then dispatch:
prop_resamplemap!(..., ::Python334Policy) = ...
prop_resamplemap!(..., ::CorrectedPolicy) = ...
```

Notes:
- Keep `compat_mode::Symbol` only in constructor arguments; do not store or propagate raw symbol in runtime paths.
- After resolution, pass `policy::CompatPolicy` (or embed in context struct) through call chains.

## 3. Backend Registration Matrix
- CPU `Array`:
  - FFT: `FFTW.jl`
  - Interp/resample: reference CPU kernels (plus accelerated kernels where available)
- GPU `CuArray`:
  - FFT: `CUDA.CUFFT`
  - Interp/resample: `KernelAbstractions.jl` / `AcceleratedKernels.jl` kernels where implemented
- Future backends:
  - oneAPI/ROCm backends may be added without API break if they satisfy trait contracts

## 4. Error Handling
- Missing backend method behavior:
  - `MethodError`
- Unsupported feature behavior:
  - explicit `ArgumentError`/`ErrorException` with supported alternatives

## 5. Contract Tests
- [ ] Trait dispatch tests by backend
- [ ] FFT equivalence tests across backends
- [ ] Interpolation consistency tests across backends
- [ ] Explicit no-scalar-indexing checks on GPU paths
- [ ] Inference checks for hot kernels (`@code_warntype`/equivalent CI gate)
