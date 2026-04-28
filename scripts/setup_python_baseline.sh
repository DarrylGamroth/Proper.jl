#!/usr/bin/env bash
set -euo pipefail

url="https://sourceforge.net/projects/proper-library/files/proper_v3.3.4_python.zip/download"
target="${PYPROPER_ROOT:-${1:-../proper_v3.3.4_python}}"

if [[ -d "${target}/proper" ]]; then
  echo "Python PROPER baseline already available: ${target}"
  exit 0
fi

tmpdir="$(mktemp -d)"
trap 'rm -rf "${tmpdir}"' EXIT

mkdir -p "$(dirname "${target}")"
curl -L --retry 3 -o "${tmpdir}/proper_v3.3.4_python.zip" "${url}"
unzip -q "${tmpdir}/proper_v3.3.4_python.zip" -d "${tmpdir}"

rm -rf "${target}"
mv "${tmpdir}/proper_v3.3.4_python" "${target}"
test -d "${target}/proper"

echo "Python PROPER baseline installed: ${target}"
