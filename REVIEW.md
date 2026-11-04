# PROPER 0.4.2 accumulated propagation review

Review candidate: `codex/prepared-propagation`, implementation `1d4d16e`.
Release base: `v0.4.1` (`137abb2`). Cleanup `654c696` was fast-forwarded onto
performance head `789374d`; its paired-kernel changes preserve parsed expressions
and operand order. The candidate version is unpublished **0.4.2**.

## Scope and preserved contract

The new public `PreparedPropagation`, `prepare_propagation` and `propagate!`
provide fixed-geometry, envelope-only PTP/WTS/STW field transforms. Preparation
retains phase data and shared FFT scratch/plans; execution does not advance
wavefront geometry. Public input/output fields remain separate from FFT scratch.
Callers serialize scratch users. Dense and strided output storage must have
independent writable elements. Read-only sampled multipliers may use indexed
views. The ordinary `prop_*` routines remain the numerical reference.

Compact quadratic caches copy already evaluated factors through exact reflection.
The separately qualified CUDA spatial chirp policy recomputes the original
expression. Private raw-order, two-field and terminal-transfer kernels preserve
per-field FFT normalization and multiplication order. Private prepared phase
arithmetic retains FP32 fields and a bounded FP32 residual/FMA path, with the
original Float64 phase and promoted field-product fallback outside its domain.
No fast-math, FTZ, optical coefficient, coordinate or calibration change is made.

Public CPU/CUDA/AMDGPU operation remains portable; CUDA performance policies
have conservative fallbacks. New hardware needs numerical and measured performance
qualification. Julia compatibility remains **1.10**. SparseArrays was added only
to test extras for the output-storage regression, not runtime dependencies.

## Independent review and disposition

A separate read-only Sol review assessed the entire `v0.4.1` accumulated diff,
including propagation equations, compact/odd-grid indexing, phase fallback,
visitor behavior, shared storage and coupled AOS use.

**PRP-OUT-01 — Medium, high confidence, resolved.** Shape/backend/precision
checks admitted sparse output matrices, whose parallel insertion can mutate
shared CSC indexing storage. Static reproduction confirmed admission without
launching unsafe work. The public boundary now admits dense/strided storage and
rejects zero strides on non-singleton axes. Sparse, indexed and zero-stride
outputs fail before pack; poisoned FFT scratch detects premature execution.
Dense and positive/negative strided outputs are compared, including capture.
The reviewer checked the final fix and tests, including CUDA/AMD type admission.
No remaining source-review blocker; execution results are recorded below.

## Verification

Fresh local gates on implementation `1d4d16e`, paired with AOS `e931beb6`:

| Gate | Result |
| --- | --- |
| Full `Pkg.test`, Julia 1.10.12, two threads | Passed, including Aqua; 3,640 testset assertions. |
| Full `Pkg.test`, Julia 1.13.3, two threads | Passed, including Aqua; 3,654 testset assertions. |
| Captured CUDA propagation / prepared phase | Passed; 1,582 / 705 testset assertions. |
| Captured AMD propagation / prepared phase | Passed; 1,587 / 705 testset assertions; CUDA remained unloaded. |
| AOS compiled-path / train CPU gates, Julia 1.12.7 | Passed; 4,431 / 92 testset assertions. |
| AOS CUDA compiled paths | Passed; 6,659 testset assertions. |
| Frozen fitted SCC / LLOWFS captured replay | Passed; 54 / 27 assertions, three changing states and reset/replay. |

CPU cohorts used KA 0.10.4 / AK 0.5.11. CUDA used RTX 3080, CUDA 6.4.2 /
KA 0.9.45 / AK 0.5.11. GPU execution used captured graphs, one Julia/BLAS
thread and serial jobs. CPU assertion-count differences are in the unchanged
version-dependent run-entrypoint tests. Sparse/zero-stride admission tests and
valid strided outputs passed on both CPU runtimes and CUDA.

Reproduce with a temporary environment developing the candidate source (add
CUDA or AMDGPU plus FFTW, KernelAbstractions and Test for direct GPU selectors):

```sh
JULIA_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1 julia +1.10 --startup-file=no --project=ENV -e 'using Pkg; Pkg.test("Proper"; julia_args=["--threads=2"])'
JULIA_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1 julia --startup-file=no --project=ENV -e 'using Pkg; Pkg.test("Proper"; julia_args=["--threads=2"])'
JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 julia --startup-file=no -O2 --project=ENV test/test_prepared_propagation.jl --cuda
JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 julia --startup-file=no -O2 --project=ENV test/test_prepared_phase.jl --cuda
```

For the local AMD cohort, use AMDGPU 2.9.2 / KA 0.10.1 / AK 0.5.11 and
`HSA_OVERRIDE_GFX_VERSION=10.3.0`; select `AMDGPU.HIPDevice(1)` before including
the same selectors with `--amdgpu`. That override is local qualification setup,
not a package policy. Companion AOS `REVIEW.md` records its fitted and exposure
commands, errors and frozen asset identities. Temporary source environments/logs
are removed after recording results; no overrides enter application manifests.

The AOS frozen-model checks use independent original numerical executors and
unchanged tolerances. Prior cleanup validation also covered AMD prepared fields,
phase and compiled/trained AOS callers. No new Python upstream parity corpus or
Roman/WFIRST exhaustive matrix run is claimed here. Julia-minimum CPU checks
qualify that runtime surface only, not minimum-version GPU library combinations.

## Coupled release and limits

Release this **0.4.2** candidate before companion AOS **0.3.20**. The AOS branch
uses the new public API plus private train/batch/cancellation lowerings, so review
and validate those packages together. Inspect the whole change with
`git diff v0.4.1..codex/prepared-propagation`; do not review cleanup alone.
Maintained API, numerical and prepared-execution guides describe the contracts.

Application model/calibration assets, source pins and deployed services are
untouched. Performance evidence belongs to the instrument investigation:
latest bounded matched source result **43.532471 updates/s per camera** and
**22.298623 ms** median, with 50 Hz still unmet. No new throughput or latency
qualification is claimed. No tags, publication or registration were performed;
human release approval and any additional platform/external-reference release
gates remain pending.
