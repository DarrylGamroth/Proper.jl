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
  echo "Python setup interpreter not found: ${PYTHON_BIN}" >&2
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

show_help() {
  cat <<'EOF'
Usage: ./scripts/profile_wfirst_phaseb_cpu.sh [--cases a,b,c] [--samples N]

Options:
  --cases     Comma-separated WFIRST case names
  --samples   Timing samples per case

Environment overrides:
  WFIRST_CASES
  WFIRST_SAMPLES
EOF
}

cases_csv="${WFIRST_CASES:-compact_hlc,full_hlc,compact_spc_spec_long,full_spc_spec_long}"
samples="${WFIRST_SAMPLES:-7}"

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
echo "Using WFIRST_CASES=${cases_csv}"
echo "Using WFIRST_SAMPLES=${samples}"

IFS=',' read -r -a cases <<<"${cases_csv}"
for case_name in "${cases[@]}"; do
  echo "[profile] Julia WFIRST Phase B ${case_name//_/ }"
  bench_julia bench/julia/wfirst_phaseb/run_case.jl \
    --case "${case_name}" \
    --data-root "${WFIRST_PHASEB_DATA_ROOT}" \
    --samples "${samples}" >/dev/null
done

bench_julia bench/julia/wfirst_phaseb/profile_cpu_cases.jl --cases "${cases_csv}"
