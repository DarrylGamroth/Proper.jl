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
Python parity environment is not ready for the WFIRST Phase B CPU benchmark.
Selected interpreter: ${PYTHON_BIN}

Recommended fix:
  ./scripts/setup_parity_venv.sh
EOF
  exit 1
fi

if [[ -z "${JULIA_NUM_THREADS:-}" ]]; then
  JULIA_NUM_THREADS="$(python3 - <<'PY'
import os
print(max(2, min(8, os.cpu_count() or 2)))
PY
)"
  export JULIA_NUM_THREADS
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

setup_json="$(${PYTHON_BIN} bench/python/setup_wfirst_phaseb_compat_data.py)"
export WFIRST_PHASEB_DATA_ROOT="$(SETUP_JSON="${setup_json}" python3 - <<'PY'
import json, os
print(json.loads(os.environ['SETUP_JSON'])['compat_root'])
PY
)"

echo "Using WFIRST_PHASEB_DATA_ROOT=${WFIRST_PHASEB_DATA_ROOT}"
echo "Using JULIA_NUM_THREADS=${JULIA_NUM_THREADS}"

cases=(
  compact_hlc
  full_hlc
  compact_spc_spec_long
  full_spc_spec_long
  compact_spc_wide
  full_spc_wide
)

for case_name in "${cases[@]}"; do
  pretty_case="${case_name//_/ }"
  run_step "Python WFIRST Phase B ${pretty_case}" "${PYTHON_BIN}" bench/python/wfirst_phaseb_external.py --case "${case_name}" --write-output-prefix "bench/reports/python_wfirst_phaseb_${case_name}"
  run_step "Julia WFIRST Phase B ${pretty_case}" julia --project=. bench/julia/wfirst_phaseb/run_case.jl --case "${case_name}" --data-root "${WFIRST_PHASEB_DATA_ROOT}"
done

julia --project=. bench/julia/wfirst_phaseb/compare_cpu.jl
