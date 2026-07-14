#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PROPER_REPO_ROOT="$(cd -- "${SCRIPT_DIR}/.." && pwd)"
backend="${1:-cpu}"

if [[ "${backend}" == "all" ]]; then
  "${BASH_SOURCE[0]}" cpu
  "${BASH_SOURCE[0]}" amdgpu
  "${BASH_SOURCE[0]}" cuda
  exit 0
fi

case "${backend}" in
  cpu)
    export PROPER_BENCH_PROJECT="${PROPER_BENCH_PROJECT:-${PROPER_REPO_ROOT}/bench}"
    ;;
  amdgpu)
    export PROPER_BENCH_PROJECT="${PROPER_BENCH_PROJECT:-${PROPER_REPO_ROOT}/bench/amdgpu}"
    ;;
  cuda)
    export PROPER_BENCH_PROJECT="${PROPER_BENCH_PROJECT:-${PROPER_REPO_ROOT}/bench/cuda}"
    ;;
  *)
    echo "usage: $0 [cpu|amdgpu|cuda|all]" >&2
    exit 2
    ;;
esac

source "${SCRIPT_DIR}/bench_env.sh"

run_latency() {
  local label="$1"
  shift
  echo "[latency] ${label}"
  if [[ "${BENCH_VERBOSE:-0}" == "1" ]]; then
    "$@"
  else
    "$@" >/dev/null
  fi
}

if [[ "${backend}" == "cpu" ]]; then
  run_latency "CPU prepared-call distribution" \
    latency_julia bench/julia/latency/cpu.jl
  exit 0
fi

ensure_bench_env
probe_output=""
if probe_output=$(bench_julia "bench/julia/${backend}/probe.jl" 2>&1); then
  run_latency "${backend^^} prepared-call distribution" \
    latency_julia "bench/julia/latency/${backend}.jl"
else
  echo "[latency] ${backend^^} unavailable; writing a skipped report"
  PROPER_LATENCY_BACKEND="${backend}" \
  PROPER_LATENCY_SKIP_REASON="${probe_output}" \
    bench_julia bench/julia/latency/write_skipped.jl >/dev/null
fi
