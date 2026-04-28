#!/usr/bin/env bash
set -euo pipefail

repo="${WFIRST_MODELS_REPO:-https://github.com/ajeldorado/proper-models.git}"
ref="${WFIRST_MODELS_REF:-master}"
target="${WFIRST_MODELS_ROOT:-${1:-../proper-models}}"
python_pkg="${target}/wfirst_cgi/models_phaseb/python/wfirst_phaseb_proper"
fetch_lfs="${WFIRST_MODELS_LFS:-0}"

if [[ -d "${python_pkg}" ]]; then
  echo "WFIRST proper-models baseline already available: ${target}"
  exit 0
fi

if ! command -v git >/dev/null 2>&1; then
  echo "git is required to fetch the WFIRST proper-models baseline" >&2
  exit 1
fi

mkdir -p "$(dirname "${target}")"
rm -rf "${target}"
mkdir -p "${target}"

git -C "${target}" init -q
git -C "${target}" remote add origin "${repo}"
if [[ "${fetch_lfs}" == "1" ]]; then
  if ! command -v git-lfs >/dev/null 2>&1; then
    echo "WFIRST_MODELS_LFS=1 requires git-lfs" >&2
    exit 1
  fi
  git -C "${target}" lfs install --local
  git -C "${target}" fetch --depth 1 origin "${ref}"
  git -C "${target}" checkout -q --detach FETCH_HEAD
  git -C "${target}" lfs pull
else
  git_lfs_disabled=(
    -c filter.lfs.process=
    -c filter.lfs.smudge=
    -c filter.lfs.clean=
    -c filter.lfs.required=false
  )
  GIT_LFS_SKIP_SMUDGE=1 git -C "${target}" "${git_lfs_disabled[@]}" fetch --depth 1 origin "${ref}"
  GIT_LFS_SKIP_SMUDGE=1 git -C "${target}" "${git_lfs_disabled[@]}" checkout -q --detach FETCH_HEAD
fi

test -d "${python_pkg}"

echo "WFIRST proper-models baseline installed: ${target}"
git -C "${target}" rev-parse HEAD
