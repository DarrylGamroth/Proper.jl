#!/usr/bin/env bash
set -euo pipefail

source "$(dirname "$0")/benchmark_cuda_lib.sh"
source "$(dirname "$0")/benchmark_amdgpu_lib.sh"

run_step() {
  local label="$1"
  shift

  echo "[bench] ${label}"
  if [[ "${BENCH_VERBOSE:-0}" == "1" ]]; then
    "$@"
  else
    "$@" >/dev/null
  fi
}

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
  ./scripts/setup_python_baseline.sh
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
run_step "Python steady-state baseline" "${PYTHON_BIN}" bench/python/run.py
run_step "Julia steady-state workload" julia --project=. bench/julia/steady_state/run.jl
run_step "Julia supported CPU kernels" julia --project=. bench/julia/steady_state/supported_kernels.jl
run_step "Julia phase-2 kernels" julia --project=. bench/julia/steady_state/phase2_kernels.jl
run_step "Julia refactor kernels" julia --project=. bench/julia/steady_state/refactor_kernels.jl
run_step "Julia KA interpolation pilot" julia --project=. bench/julia/steady_state/ka_interp_kernels.jl
run_step "Julia KA geometry/sampling pilot" julia --project=. bench/julia/steady_state/ka_geometry_sampling_kernels.jl
run_step "Julia example workflows" julia --project=. bench/julia/steady_state/example_workflows.jl
if [[ "${BENCH_INCLUDE_WFIRST_CPU:-0}" == "1" ]]; then
  run_step "WFIRST Phase B CPU comparison" ./scripts/benchmark_wfirst_phaseb_cpu.sh
fi
run_cuda_benchmarks
run_amdgpu_benchmarks
run_step "Julia cold-start / TTFx" julia --project=. bench/julia/cold_start/run.jl
julia --project=. bench/reports/summarize.jl
