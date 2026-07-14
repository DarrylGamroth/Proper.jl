#!/usr/bin/env bash
set -euo pipefail

source "$(dirname "$0")/bench_env.sh"

ensure_bench_env

workload="${PROPER_THREAD_TOPOLOGY_WORKLOAD:-core}"
case "${workload}" in
  core)
    default_matrix="1:1 1:2 1:4 1:8 2:4 4:2 8:1"
    ;;
  batch)
    default_matrix="1:1 2:1 4:1 8:1 1:4 1:8 2:4 4:2"
    ;;
  *)
    echo "invalid workload '${workload}'; expected core or batch" >&2
    exit 2
    ;;
esac

matrix="${PROPER_THREAD_TOPOLOGY_MATRIX:-${default_matrix}}"
final_report="${PROPER_BENCH_REPORT:-${PROPER_REPO_ROOT}/bench/reports/julia_thread_topology_cpu_${workload}.json}"
scratch_dir="$(mktemp -d "${TMPDIR:-/tmp}/proper-thread-topology.XXXXXX")"
reference_path="${scratch_dir}/reference.jls"
case_reports=()
read -r -a topologies <<< "${matrix}"

if [[ ${#topologies[@]} -eq 0 ]]; then
  echo "thread topology matrix must not be empty" >&2
  exit 2
fi

cleanup() {
  rm -rf "${scratch_dir}"
}
trap cleanup EXIT

case_index=0
for topology in "${topologies[@]}"; do
  if [[ ! "${topology}" =~ ^([1-9][0-9]*):([1-9][0-9]*)$ ]]; then
    echo "invalid thread topology '${topology}'; expected JULIA_THREADS:FFTW_THREADS" >&2
    exit 2
  fi

  julia_threads="${BASH_REMATCH[1]}"
  fftw_threads="${BASH_REMATCH[2]}"
  case_index=$((case_index + 1))
  case_report="${scratch_dir}/case${case_index}_j${julia_threads}_f${fftw_threads}.json"
  case_reports+=("${case_report}")
  echo "[bench] Julia threads=${julia_threads}, FFTW threads=${fftw_threads}, BLAS threads=${PROPER_BENCH_BLAS_THREADS:-1}"

  (
    export PROPER_BENCH_JULIA_THREADS="${julia_threads}"
    export PROPER_BENCH_FFTW_THREADS="${fftw_threads}"
    export PROPER_BENCH_BLAS_THREADS="${PROPER_BENCH_BLAS_THREADS:-1}"
    export PROPER_BENCH_WORKLOAD="${workload}"
    export PROPER_BENCH_REFERENCE="${reference_path}"
    export PROPER_BENCH_REPORT="${case_report}"
    bench_julia --threads="${julia_threads}" "${PROPER_REPO_ROOT}/bench/julia/steady_state/thread_topology_case.jl"
  )
done

(
  export PROPER_BENCH_REPORT="${final_report}"
  bench_julia "${PROPER_REPO_ROOT}/bench/julia/steady_state/aggregate_thread_topology.jl" "${case_reports[@]}"
)

echo "[bench] wrote ${final_report}"
