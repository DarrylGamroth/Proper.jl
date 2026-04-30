#!/usr/bin/env bash
set -euo pipefail

source "$(dirname "$0")/bench_env.sh"
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
  ./scripts/setup_python_baseline.sh
  ./scripts/setup_wfirst_models_baseline.sh
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

show_help() {
  cat <<'EOF'
Usage: ./scripts/benchmark_wfirst_phaseb_cpu.sh [--cases a,b,c] [--samples N] [--parity-only]

Options:
  --cases         Comma-separated WFIRST case names to run
  --samples       Timing samples per case in timed mode
  --parity-only   Skip timing loops and run correctness/output generation only

Environment overrides:
  WFIRST_CASES
  WFIRST_SAMPLES
  WFIRST_PARITY_ONLY=1
EOF
}

cases_csv="${WFIRST_CASES:-}"
samples="${WFIRST_SAMPLES:-3}"
parity_only="${WFIRST_PARITY_ONLY:-0}"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --cases)
      cases_csv="$2"
      shift 2
      ;;
    --samples)
      samples="$2"
      shift 2
      ;;
    --parity-only)
      parity_only=1
      shift
      ;;
    --help|-h)
      show_help
      exit 0
      ;;
    *)
      echo "Unknown argument: $1" >&2
      show_help >&2
      exit 1
      ;;
  esac
done

setup_json="$(${PYTHON_BIN} bench/python/setup_wfirst_phaseb_compat_data.py)"
export WFIRST_PHASEB_DATA_ROOT="$(SETUP_JSON="${setup_json}" python3 - <<'PY'
import json, os
print(json.loads(os.environ['SETUP_JSON'])['compat_root'])
PY
)"

echo "Using WFIRST_PHASEB_DATA_ROOT=${WFIRST_PHASEB_DATA_ROOT}"
echo "Using JULIA_NUM_THREADS=${JULIA_NUM_THREADS}"

if [[ -z "${cases_csv}" ]]; then
  cases=(
    compact_hlc
    full_hlc
    compact_spc_spec_short
    full_spc_spec_short
    compact_spc_ifs_short
    full_spc_ifs_short
    compact_spc_spec_long
    full_spc_spec_long
    compact_spc_ifs_long
    full_spc_ifs_long
    compact_spc_wide
    full_spc_wide
  )
else
  IFS=',' read -r -a cases <<<"${cases_csv}"
fi

python_extra_args=(--samples "${samples}")
julia_extra_args=(--samples "${samples}")
if [[ "${parity_only}" == "1" ]]; then
  python_extra_args+=(--parity-only)
  julia_extra_args+=(--parity-only)
fi

for case_name in "${cases[@]}"; do
  pretty_case="${case_name//_/ }"
  run_step "Python WFIRST Phase B ${pretty_case}" "${PYTHON_BIN}" bench/python/wfirst_phaseb_external.py --case "${case_name}" --write-output-prefix "bench/reports/python_wfirst_phaseb_${case_name}" "${python_extra_args[@]}"
  run_step "Julia WFIRST Phase B ${pretty_case}" bench_julia bench/julia/wfirst_phaseb/run_case.jl --case "${case_name}" --data-root "${WFIRST_PHASEB_DATA_ROOT}" "${julia_extra_args[@]}"
done

"${PYTHON_BIN}" bench/python/compare_wfirst_phaseb_outputs.py --cases "$(IFS=,; echo "${cases[*]}")"
