#!/usr/bin/env bash
set -euo pipefail

python3 bench/python/run.py
julia --project=. bench/julia/steady_state/run.jl
julia --project=. bench/julia/cold_start/run.jl
julia --project=. bench/reports/summarize.jl
