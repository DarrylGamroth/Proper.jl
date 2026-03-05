#!/usr/bin/env bash
set -euo pipefail

python3 bench/python/run.py
julia --project=. bench/julia/steady_state/run.jl
julia --project=. bench/julia/steady_state/phase2_kernels.jl
julia --project=. bench/julia/steady_state/refactor_kernels.jl
julia --project=. bench/julia/steady_state/example_workflows.jl
julia --project=. bench/julia/cold_start/run.jl
julia --project=. bench/reports/summarize.jl
