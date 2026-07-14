#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
python_bin="${PYTHON_BIN:-${repo_root}/.venv-parity/bin/python}"

if [[ ! -x "${python_bin}" ]]; then
  echo "Python parity interpreter not found: ${python_bin}" >&2
  exit 1
fi

scratch_dir="$(mktemp -d "${TMPDIR:-/tmp}/proper-parity-repro.XXXXXX")"
trap 'rm -rf "${scratch_dir}"' EXIT

for run in first second; do
  output_dir="${scratch_dir}/${run}"
  mkdir -p "${output_dir}"
  SOURCE_DATE_EPOCH=0 "${python_bin}" \
    "${repo_root}/test/parity/generate_python_baseline.py" \
    --output-dir "${output_dir}" >/dev/null
  SOURCE_DATE_EPOCH=0 "${python_bin}" \
    "${repo_root}/test/parity/generate_python_example_metrics.py" \
    --output-dir "${output_dir}" >/dev/null
done

diff -ru "${scratch_dir}/first" "${scratch_dir}/second"
echo "Python parity baseline generation is byte-reproducible."
