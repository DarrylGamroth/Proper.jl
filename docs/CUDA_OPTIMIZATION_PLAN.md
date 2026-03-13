# CUDA Optimization Plan

Date: 2026-03-13  
Owner: Proper.jl port effort  
Scope: steady-state CUDA runtime for parity-preserving workloads and supported GPU kernels

## Problem Statement
- CUDA is functional, but current steady-state GPU performance is inconsistent.
- Interpolation kernels are already competitive on GPU.
- Propagation and some geometry/mask kernels are still limited by:
  - host-staged temporary map generation
  - CPU-owned workspace buffers reused from GPU paths
  - extra FFT temporaries on the CUDA path
  - unconditional kernel synchronization in low-level helpers
  - full-frame launches for small geometric masks
  - mixed FP32/FP64 benchmark interpretation

## Current External Baseline
Reference benchmark: user run on 2026-03-13, RTX 3050 Ti Laptop GPU.

- Steady-state workload:
  - Python CPU: `96.04 ms`
  - Julia CPU: `24.93 ms`
  - Julia CUDA: `30.52 ms`
- Supported kernel highlights:
  - GPU wins:
    - `prop_magnify_quick_mutating`: `5.68x`
    - `prop_resamplemap_mutating`: `12.14x`
    - `prop_rotate_mutating`: `1.17x`
    - `prop_szoom_mutating`: `1.32x`
  - GPU losses:
    - `prop_ptp`: `0.55x`
    - `prop_end_mutating`: `0.54x`
    - `prop_rectangle_mutating`: `0.08x`
    - `prop_rounded_rectangle_mutating`: `0.07x`

## Design Rules
- Preserve parity behavior against the patched Python 3.3.4 baseline.
- Prefer typed dispatch and backend traits over runtime branching.
- Keep CPU fast paths intact unless the GPU-oriented implementation is neutral or better on CPU.
- Avoid host staging on CUDA paths unless no equivalent device-safe algorithm exists.
- Synchronize at API/benchmark boundaries, not inside low-level kernels.

## Phases

### C0: Baseline and Bottleneck Classification
Status: Completed
- Capture CUDA steady-state and kernel-level comparison reports.
- Classify regressions by:
  - host staging / copies
  - scratch/FFT allocation
  - launch/synchronization overhead
  - precision regime differences

Outcome:
- Main issues are allocation/copy/synchronization related, not type instability.

### C1: Remove Host-Staged Phase Maps and Internal Sync Overhead
Status: Completed
- Replace broadcast/map construction on CUDA in:
  - `prop_qphase`
  - `prop_ptp`
- Move quadratic phase and frequency-domain phase application to direct KA kernels.
- Remove unconditional `AK.synchronize` from low-level KA helpers.
- Restrict rectangle/rounded-rectangle CUDA kernels to bounding-box launch regions.

Acceptance:
- No scalar indexing on CUDA paths.
- CPU tests stay green.
- CPU benchmark suite remains green.

### C2: Backend-Aware Workspace State
Status: In Progress
- Make workspace storage backend-aware instead of CPU-only.
- Preserve device buffers for:
  - interpolation axes
  - mask scratch
  - FFT scratch / cached maps
- Ensure `WaveFront` and `RunContext` can carry backend-appropriate scratch state.

Completed in this slice:
- `InterpWorkspace` is now backend-aware.
- `MaskWorkspace.mask` is now backend-aware.
- `WaveFront` and `RunContext` now construct workspaces from the field/output backend type.
- `prop_magnify!` and `prop_resamplemap!` now fill interpolation axes through typed loop/KA dispatch instead of host scalar loops.

Remaining in C2:
- backend-aware FFT/cache state remains deferred to `C3`
- polygon vertex scratch remains host-owned because the vectors are small and not a measured hotspot

Targets:
- Eliminate repeated `backend_adapt(...)` staging on hot CUDA paths.
- Provide a clean basis for device-side cache reuse.

### C3: CUDA FFT Scratch and Propagation Reuse
Status: In Progress
- Add CUDA-oriented scratch reuse for:
  - `prop_ptp`
  - `prop_wts`
  - `prop_stw`
  - `prop_end`
- Keep dispatch explicit by backend/FFT style.
- Minimize temporary FFT arrays and repeated device allocations per call.

Completed in this slice:
- backend-aware complex FFT scratch is now preserved in `FFTWorkspace`
- CUDA propagation paths for `prop_ptp`, `prop_wts`, and `prop_stw` now use workspace-backed in-place `fft!` / `bfft!` flows
- CPU planned FFT path is preserved unchanged behind typed dispatch

Remaining in C3:
- validate CUDA in-place FFT behavior and timings on hardware
- decide whether `prop_end` needs additional backend scratch/state or whether current KA copy path is sufficient

Targets:
- Improve steady-state propagation on GPU.
- Reduce gap between CPU and CUDA in the end-to-end workload.

### C4: Remove Remaining Host Staging in Mask/Map Paths
Status: Pending
- Audit and convert remaining CUDA-visible host staging in:
  - `prop_circular_aperture`
  - other mask wrappers using CPU workspace scratch
  - map-generation helpers that still materialize host intermediates

Targets:
- Avoid device-host-device round-trips in common aperture workflows.

### C5: Benchmark Clarity and Precision Split
Status: Pending
- Split CUDA benchmark reporting by precision regime where useful:
  - wavefront FP64 propagation
  - image/interpolation FP32 kernels
- Add clearer comparison tables for:
  - CPU vs CUDA by kernel family
  - FP32 vs FP64 where behavior differs materially

Targets:
- Make GPU results interpretable without conflating precision and algorithm effects.

## Execution Log
- 2026-03-13: Plan created from CUDA benchmark audit and RTX 3050 Ti report.
- 2026-03-13: C1 implemented.
  - Added direct CUDA KA kernels for quadratic phase and frequency-domain phase application.
  - Routed CUDA `prop_qphase` and `prop_ptp` through typed execution selectors.
  - Removed unconditional low-level KA synchronization.
  - Restricted rectangle and rounded-rectangle CUDA launches to computed bounding boxes.
  - Added CUDA smoke coverage for `prop_ptp`.
- 2026-03-13: C2 slice 1 implemented.
  - `InterpWorkspace` now preserves the requested backend.
  - `MaskWorkspace.mask` now preserves the wavefront/backend type.
  - `WaveFront` and `RunContext` now build backend-aware workspace state for interpolation and mask paths.
  - added device-safe affine-axis fill for interpolation coordinates and CUDA smoke coverage for `prop_resamplemap!`.
  - local validation:
    - `julia --project=. test/runtests.jl`: pass
    - `./scripts/benchmark_all.sh`: pass on non-CUDA machine, CUDA lane skipped cleanly
  - local benchmark snapshot:
    - `simple_prescription_256`: `4.01 MiB -> 3.50 MiB`
    - `simple_telescope_256`: `5.02 MiB -> 3.50 MiB`
    - `psdtest_128`: `1.64 MiB -> 1.51 MiB`
- 2026-03-13: C3 slice 1 implemented.
  - `FFTWorkspace` now preserves backend complex scratch buffers while keeping CPU FFT plan caching intact.
  - CUDA propagation routes for `prop_ptp`, `prop_wts`, and `prop_stw` now use workspace-backed in-place `fft!` / `bfft!` flows.
  - added CPU scratch-reuse regression coverage and CUDA smoke coverage for backend-preserving FFT scratch.
  - local validation:
    - `julia --project=. test/runtests.jl`: pass
    - `./scripts/benchmark_all.sh`: pass on non-CUDA machine, CUDA lane skipped cleanly
  - local benchmark snapshot:
    - steady-state CPU: `27.89 ms -> 26.92 ms`
    - phase-2 `prop_ptp`: `12.85 ms -> 10.76 ms`
    - example allocations unchanged from C2 slice (`3.50 MiB` / `3.50 MiB` / `1.51 MiB`)
- 2026-03-13: Awaiting rerun on CUDA hardware for updated C1/C2/C3 benchmark data.
