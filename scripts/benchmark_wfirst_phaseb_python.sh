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

if ! command -v "${PYTHON_BIN}" >/dev/null 2>&1; then
  echo "Python benchmark interpreter not found: ${PYTHON_BIN}" >&2
  exit 1
fi

required_python_modules() {
  "$1" - <<'PY'
import astropy
import matplotlib
import numpy
import scipy
PY
}

if ! required_python_modules "${PYTHON_BIN}" >/dev/null 2>&1; then
  cat >&2 <<EOF
Python parity environment is not ready for the external WFIRST Phase B benchmark.
Selected interpreter: ${PYTHON_BIN}

Recommended fix:
  ./scripts/setup_python_baseline.sh
  ./scripts/setup_wfirst_models_baseline.sh
  ./scripts/setup_parity_venv.sh
EOF
  exit 1
fi

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

run_step "Python WFIRST Phase B compact HLC" "${PYTHON_BIN}" bench/python/wfirst_phaseb_external.py --case compact_hlc
run_step "Python WFIRST Phase B full HLC" "${PYTHON_BIN}" bench/python/wfirst_phaseb_external.py --case full_hlc
julia --project=. bench/reports/summarize.jl
