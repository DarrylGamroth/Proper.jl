#!/usr/bin/env bash
set -euo pipefail

PYTHON_BIN="${PYTHON_BIN:-}"
if [[ -z "${PYTHON_BIN}" ]]; then
  if [[ -x ".venv-parity/bin/python" ]]; then
    PYTHON_BIN=".venv-parity/bin/python"
  else
    PYTHON_BIN="python3"
  fi
fi

"${PYTHON_BIN}" bench/python/run.py
julia --project=. bench/julia/steady_state/run.jl
julia --project=. bench/julia/steady_state/phase2_kernels.jl
julia --project=. bench/julia/steady_state/refactor_kernels.jl
julia --project=. bench/julia/steady_state/example_workflows.jl
julia --project=. bench/julia/cold_start/run.jl
julia --project=. bench/reports/summarize.jl
