# Backend Traits Contract

## Purpose
Define stable trait-dispatched interfaces for array backends and computational kernels.

## Decision Links
- `D-0008` GPU scope for phase 1 (accepted)
- `D-0009` backend traits and interfaces (accepted)
- `D-0012` `prop_run_multi` semantics (accepted)
- `D-0035` remove runtime compatibility modes (accepted)
- `D-0023` type-stable core structs (accepted)
- `D-0028` dynamic dispatch minimization (accepted)

## Status
- [x] Draft pre-filled from proposed defaults
- [x] Trait interfaces frozen
- [x] Implementations registered for current CPU baseline and context-routed interpolation/propagation paths
- [x] CUDA package extension registers `CuArray` backend, FFT, and KA kernel traits without making `CUDA.jl` a hard dependency
- [x] AMDGPU package extension registers `ROCArray` backend, rocFFT, and KA kernel traits without making `AMDGPU.jl` a hard dependency

## 1. Trait Model
Proposed trait families:
- `BackendStyle(::Type{<:AbstractArray})`
- `FFTStyle(::Type{<:AbstractArray})`
- `InterpStyle(::Type{<:AbstractArray})`
- `RNGStyle(::Type)`

Compatibility requirements:
- Traits are dispatch aids; avoid runtime mutable global backend switches.
- Methods should preserve array backend where practical.
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

### 2.5 Baseline Contract
- Parity behavior is anchored to the patched Python 3.3.4 executable baseline.
- Do not add runtime compatibility mode flags in context or kernel call chains.
- Divergences from baseline must be documented and test-covered in parity reports.

## 3. Backend Registration Matrix
- CPU `Array`:
  - FFT: `FFTW.jl`
  - Interp/resample: reference CPU kernels (plus accelerated kernels where available)
- GPU `CuArray`:
  - FFT: `CUDA.CUFFT`
  - Interp/resample: `KernelAbstractions.jl` / `AcceleratedKernels.jl` kernels where implemented
  - Registration: package extension `ProperCUDAExt`
- GPU `ROCArray`:
  - FFT: `AMDGPU.rocFFT`
  - Interp/resample: `KernelAbstractions.jl` / `AcceleratedKernels.jl` kernels where implemented
  - Registration: package extension `ProperAMDGPUExt`
- Future backends:
  - oneAPI/ROCm backends may be added without API break if they satisfy trait contracts

## 4. Error Handling
- Missing backend method behavior:
  - `MethodError`
- Unsupported feature behavior:
  - explicit `ArgumentError`/`ErrorException` with supported alternatives

## 5. Contract Tests
- [x] Trait dispatch tests by backend (CPU baseline + optional CUDA smoke when CUDA is available)
- [ ] FFT equivalence tests across backends
- [x] Interpolation consistency tests across context/style-dispatched entry points
- [x] Explicit no-scalar-indexing checks on GPU smoke path (`CUDA.allowscalar(false)`, availability-gated)
- [ ] Inference checks for hot kernels (`@code_warntype`/equivalent CI gate)
