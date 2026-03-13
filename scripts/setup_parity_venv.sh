#!/usr/bin/env bash
set -euo pipefail

BOOTSTRAP_PYTHON="${BOOTSTRAP_PYTHON:-python3}"
VENV_DIR="${VENV_DIR:-.venv-parity}"
VENV_PYTHON="${VENV_DIR}/bin/python"

python_has_module() {
  "$1" - "$2" <<'PY'
import importlib.util
import sys

module_name = sys.argv[1]
sys.exit(0 if importlib.util.find_spec(module_name) is not None else 1)
PY
}

create_venv() {
  "${BOOTSTRAP_PYTHON}" -m venv "${VENV_DIR}"
}

if ! command -v "${BOOTSTRAP_PYTHON}" >/dev/null 2>&1; then
  echo "Bootstrap Python not found: ${BOOTSTRAP_PYTHON}" >&2
  exit 1
fi

if [[ ! -x "${VENV_PYTHON}" ]]; then
  create_venv
fi

if ! python_has_module "${VENV_PYTHON}" pip; then
  echo "Bootstrapping pip in ${VENV_DIR}" >&2
  "${VENV_PYTHON}" -m ensurepip --upgrade
fi

"${VENV_PYTHON}" -m pip install --upgrade pip
"${VENV_PYTHON}" -m pip install numpy scipy astropy matplotlib

echo "Parity environment ready: ${VENV_PYTHON}"
