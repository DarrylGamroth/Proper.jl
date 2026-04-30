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

show_help() {
  cat <<'EOF'
Usage: ./scripts/profile_wfirst_phaseb_stages.sh [--cases a,b,c]

Options:
  --cases     Comma-separated WFIRST case names

Environment overrides:
  WFIRST_CASES
EOF
}

cases_csv="${WFIRST_CASES:-full_spc_spec_long}"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --cases)
      cases_csv="$2"
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
echo "Using WFIRST_CASES=${cases_csv}"

bench_julia bench/julia/wfirst_phaseb/profile_stage_cases.jl \
  --cases "${cases_csv}" \
  --data-root "${WFIRST_PHASEB_DATA_ROOT}"
