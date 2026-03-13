#!/usr/bin/env bash
set -euo pipefail

required_python_modules() {
  "$1" - <<'PY'
import astropy
import matplotlib
import numpy
import scipy
PY
}

print_python_setup_help() {
  cat >&2 <<EOF
Python parity environment is not ready for benchmarks.
Selected interpreter: ${PYTHON_BIN}

Required packages:
  numpy scipy astropy matplotlib

Recommended fix:
  ./scripts/setup_parity_venv.sh

Alternative:
  ${PYTHON_BIN} -m pip install --upgrade pip numpy scipy astropy matplotlib

Interpreter override:
  PYTHON_BIN=/path/to/python ./scripts/benchmark_all.sh
EOF
}

PYTHON_BIN="${PYTHON_BIN:-}"
if [[ -z "${PYTHON_BIN}" ]]; then
  if [[ -x ".venv-parity/bin/python" ]]; then
    PYTHON_BIN=".venv-parity/bin/python"
  else
    PYTHON_BIN="python3"
  fi
fi

if ! command -v "${PYTHON_BIN}" >/dev/null 2>&1; then
  echo "Python benchmark interpreter not found: ${PYTHON_BIN}" >&2
  exit 1
fi

if ! required_python_modules "${PYTHON_BIN}" >/dev/null 2>&1; then
  print_python_setup_help
  exit 1
fi

echo "Using Python benchmark interpreter: ${PYTHON_BIN}"
"${PYTHON_BIN}" bench/python/run.py
julia --project=. bench/julia/steady_state/run.jl
julia --project=. bench/julia/steady_state/phase2_kernels.jl
julia --project=. bench/julia/steady_state/refactor_kernels.jl
julia --project=. bench/julia/steady_state/ka_interp_kernels.jl
julia --project=. bench/julia/steady_state/example_workflows.jl
julia --project=. bench/julia/cuda/steady_state.jl
julia --project=. bench/julia/cuda/supported_kernels.jl
julia --project=. bench/julia/cold_start/run.jl
julia --project=. bench/reports/summarize.jl
