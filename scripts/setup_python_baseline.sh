#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
target="${PYPROPER_ROOT:-${1:-${repo_root}/../proper_v3.3.4_python}}"
seed_url="${PYPROPER_SEED_URL:-https://sourceforge.net/projects/proper-library/files/proper_v3.3.5_python.zip/download}"
seed_sha256="${PYPROPER_SEED_SHA256:-5c25bc4ca80efb088990f1d6be231fe5583a806ff806de17ddc26026f2b23d87}"
expected_snapshot="2f6f351715a49524f01aded1baedf7ef9c41bb40a3f738a4e7f481ed1daa7382"
snapshot_tool="${repo_root}/test/parity/baseline_provenance.py"
reconstruction_patch="${repo_root}/scripts/python_baseline_335_to_334.patch"

source_snapshot() {
  python3 "${snapshot_tool}" --source-snapshot "$1"
}

file_sha256() {
  if command -v sha256sum >/dev/null 2>&1; then
    sha256sum "$1" | awk '{print $1}'
  else
    shasum -a 256 "$1" | awk '{print $1}'
  fi
}

if [[ -d "${target}/proper" ]]; then
  actual_snapshot="$(source_snapshot "${target}")"
  if [[ "${actual_snapshot}" != "${expected_snapshot}" ]]; then
    echo "Python PROPER baseline at ${target} has unexpected source snapshot." >&2
    echo "Expected: ${expected_snapshot}" >&2
    echo "Actual:   ${actual_snapshot}" >&2
    echo "Move that tree aside or set PYPROPER_ROOT to the accepted patched 3.3.4 baseline." >&2
    exit 1
  fi
  echo "Verified Python PROPER 3.3.4 baseline: ${target}"
  exit 0
fi

for command in curl patch python3 unzip; do
  if ! command -v "${command}" >/dev/null 2>&1; then
    echo "Required baseline bootstrap command not found: ${command}" >&2
    exit 1
  fi
done

if [[ -e "${target}" ]]; then
  echo "Refusing to replace non-baseline path: ${target}" >&2
  exit 1
fi

tmpdir="$(mktemp -d)"
trap 'rm -rf "${tmpdir}"' EXIT
archive="${tmpdir}/proper_v3.3.5_python.zip"

curl -LfsS --retry 3 -o "${archive}" "${seed_url}"
actual_seed_sha256="$(file_sha256 "${archive}")"
if [[ "${actual_seed_sha256}" != "${seed_sha256}" ]]; then
  echo "Python PROPER bootstrap archive checksum mismatch." >&2
  echo "Expected: ${seed_sha256}" >&2
  echo "Actual:   ${actual_seed_sha256}" >&2
  exit 1
fi

unzip -q "${archive}" -d "${tmpdir}"
seed_root="${tmpdir}/proper_v3.3.5_python"
test -d "${seed_root}/proper"
patch --batch --forward -p1 -d "${seed_root}" < "${reconstruction_patch}"
rm "${seed_root}/proper/prop_set_rayleighfactor.py"
rm "${seed_root}/PROPER_manual_v3.3.5.pdf"

actual_snapshot="$(source_snapshot "${seed_root}")"
if [[ "${actual_snapshot}" != "${expected_snapshot}" ]]; then
  echo "Reconstructed Python PROPER source snapshot mismatch." >&2
  echo "Expected: ${expected_snapshot}" >&2
  echo "Actual:   ${actual_snapshot}" >&2
  exit 1
fi

mkdir -p "$(dirname "${target}")"
mv "${seed_root}" "${target}"
echo "Installed verified patched Python PROPER 3.3.4 baseline: ${target}"
