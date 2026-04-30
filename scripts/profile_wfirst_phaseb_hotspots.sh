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
Usage: ./scripts/profile_wfirst_phaseb_hotspots.sh [--case name] [--mincount N] [--maxdepth N]

Options:
  --case       WFIRST case name
  --mincount   Minimum sample count shown in the flat profile
  --maxdepth   Maximum stack depth printed

Environment overrides:
  WFIRST_CASE
EOF
}

case_name="${WFIRST_CASE:-full_spc_spec_long}"
mincount="25"
maxdepth="32"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --case)
      case_name="$2"
      shift 2
      ;;
    --mincount)
      mincount="$2"
      shift 2
      ;;
    --maxdepth)
      maxdepth="$2"
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
echo "Using WFIRST_CASE=${case_name}"

bench_julia bench/julia/wfirst_phaseb/profile_case_hotspots.jl \
  --case "${case_name}" \
  --data-root "${WFIRST_PHASEB_DATA_ROOT}" \
  --mincount "${mincount}" \
  --maxdepth "${maxdepth}"
