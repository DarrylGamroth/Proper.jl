#!/usr/bin/env bash
set -euo pipefail

BOOTSTRAP_PYTHON="${BOOTSTRAP_PYTHON:-python3}"
VENV_DIR="${VENV_DIR:-.venv-parity}"

if ! command -v "${BOOTSTRAP_PYTHON}" >/dev/null 2>&1; then
  echo "Bootstrap Python not found: ${BOOTSTRAP_PYTHON}" >&2
  exit 1
fi

if [[ ! -d "${VENV_DIR}" ]]; then
  "${BOOTSTRAP_PYTHON}" -m venv "${VENV_DIR}"
fi

"${VENV_DIR}/bin/python" -m pip install --upgrade pip
"${VENV_DIR}/bin/python" -m pip install numpy scipy astropy matplotlib

echo "Parity environment ready: ${VENV_DIR}/bin/python"
