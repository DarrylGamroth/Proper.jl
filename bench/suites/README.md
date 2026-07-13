# Benchmark Suites

This directory holds benchmark-suite helpers that are shared by focused
performance measurements. Suite helpers should stay here unless they become
general infrastructure suitable for `bench/common/`.

Benchmark entry points currently include:

- `bench/julia/steady_state/dm_projection.jl`
- `bench/julia/steady_state/zernike_synthesis.jl`
- `bench/julia/steady_state/zernike_fit.jl`
- `bench/julia/wfirst_phaseb/prepared_models.jl`

Generated reports belong under ignored `bench/reports/`.

All reports include Julia, FFTW, and BLAS thread counts. Set
`PROPER_FFTW_THREADS` before launching a benchmark when comparing FFTW thread
configurations, and keep outer Julia-thread and inner FFTW-thread experiments
separate so oversubscription is visible.
