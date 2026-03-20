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

## 4. Honest Support Matrix

This section describes the current user-facing GPU support surface. It is meant
to answer a practical question: which paths are truly backend-native today,
which require specific overloads/options, and which intentionally fail instead
of silently falling back to the host.

### 4.1 Fully Backend-Native On CUDA And AMDGPU
- Propagation core with explicit context reuse:
  - `prop_qphase`
  - `prop_ptp`
  - `prop_wts`
  - `prop_stw`
  - `prop_end!` when `out` and `wf.field` use the same backend
- GPU geometry/mask paths already exercised by the optional backend smoke:
  - `prop_circular_aperture`
  - `prop_rectangle`
  - `prop_rounded_rectangle`
- Map-application helpers when the map is already on the same backend as the
  wavefront:
  - `prop_add_phase`
  - `prop_multiply`
  - `prop_divide`
- FITS map application surface after promotion to the wavefront backend:
  - `prop_readmap`
  - `prop_errormap`
  - `prop_psd_errormap`
  - `prop_readmap` still decodes FITS on the host first, then promotes the map
    to the wavefront backend before resampling

### 4.2 Backend-Native Only For Specific Entry Points Or Options
- `prop_rotate`
  - full-image rotate paths are backend-native for the supported array
    backends
  - explicit point-sampling and coordinate-grid `prop_cubic_conv` forms are not
    the GPU performance surface
- `prop_magnify`
  - use `prop_magnify!` with an explicit `ctx` for the stable performance path
  - allocating wrappers exist for convenience but may create fresh workspace
    state
- `prop_resamplemap`
  - use `prop_resamplemap!` with an explicit `ctx` for the stable performance
    path
  - allocating wrappers exist for convenience and may allocate interpolation
    state
- `prop_end`
  - the allocating wrapper is backend-native
  - `prop_end!` is backend-native only when the output array uses the same
    backend as the wavefront field

### 4.3 Intentionally Unsupported On GPU
- `prop_cubic_conv` scalar point-sampling on GPU arrays
- `prop_cubic_conv` pointwise vector mode on GPU arrays
- `prop_cubic_conv` `grid=false` coordinate-grid mode on GPU arrays
- `prop_end!` with mismatched output and wavefront backends

These paths throw explicit `ArgumentError` rather than silently materializing
through host `Matrix(...)` fallbacks.

## 5. Error Handling
- Missing backend method behavior:
  - `MethodError`
- Unsupported feature behavior:
  - explicit `ArgumentError`/`ErrorException` with supported alternatives

## 6. Performance-Surface Contract
- For repeated GPU work, the intended performance surface is:
  - mutating `!` APIs
  - explicit `RunContext`
  - prepared execution where applicable
- Allocating wrappers remain available for convenience and portability, but they
  are not the preferred steady-state benchmark or HIL surface.
- Benchmark policy:
  - `BenchmarkTools` numbers represent warmed steady-state execution
  - cold-start / TTFx must be measured separately in fresh processes

## 7. Contract Tests
- [x] Trait dispatch tests by backend (CPU baseline + optional CUDA smoke when CUDA is available)
- [ ] FFT equivalence tests across backends
- [x] Interpolation consistency tests across context/style-dispatched entry points
- [x] Explicit no-scalar-indexing checks on GPU smoke path (`CUDA.allowscalar(false)`, availability-gated)
- [x] Warmed host-allocation regression gates for GPU propagation hot paths
- [x] Optional GPU smoke for same-backend map application, FITS map reads, error-map application, and PSD error-map application
- [ ] Inference checks for hot kernels (`@code_warntype`/equivalent CI gate)
