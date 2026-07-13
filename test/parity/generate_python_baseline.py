#!/usr/bin/env python3
import csv
import json
import os
import sys
from pathlib import Path

import numpy as np


def run_simple_case(proper):
    wf = proper.prop_begin(2.4, 0.55e-6, 256, beam_diam_fraction=0.5)
    proper.prop_circular_aperture(wf, 0.6)
    proper.prop_lens(wf, 20.0)
    proper.prop_propagate(wf, 20.0)
    psf, sampling = proper.prop_end(wf)
    return psf, float(sampling)


def complex_array_payload(array):
    """Encode a small complex matrix without requiring a binary artifact."""
    return {
        "real": np.asarray(array).real.tolist(),
        "imag": np.asarray(array).imag.tolist(),
    }


def run_wavefront_accessors_even_case(proper):
    """Exercise centered field access and wavefront addition pixel by pixel."""
    n = 8
    yy, xx = np.indices((n, n), dtype=np.float64)

    # Dyadic coefficients make the field exactly reproducible while the
    # distinct x/y terms expose transposes and incorrect center shifts.
    raw = (
        (3.0 * yy - 2.0 * xx + 1.0) / 8.0
        + 1j * (5.0 * xx - 4.0 * yy + 3.0) / 16.0
    )
    centered_addition = (
        (7.0 * yy + xx - 5.0) / 32.0
        + 1j * (2.0 * yy - 3.0 * xx + 7.0) / 8.0
    )
    scalar_addition = 0.75 - 1.25j

    wf_accessors = proper.prop_begin(1.0, 500e-9, n, beam_diam_fraction=0.5)
    wf_accessors.wfarr = raw.copy()

    wf_scalar = proper.prop_begin(1.0, 500e-9, n, beam_diam_fraction=0.5)
    wf_scalar.wfarr = raw.copy()
    proper.prop_add_wavefront(wf_scalar, scalar_addition)

    wf_matrix = proper.prop_begin(1.0, 500e-9, n, beam_diam_fraction=0.5)
    wf_matrix.wfarr = raw.copy()
    proper.prop_add_wavefront(wf_matrix, centered_addition)

    return {
        "inputs": {
            "raw_field": complex_array_payload(raw),
            "centered_addition": complex_array_payload(centered_addition),
            "scalar_addition": {
                "real": float(scalar_addition.real),
                "imag": float(scalar_addition.imag),
            },
        },
        "outputs": {
            "get_wavefront": complex_array_payload(proper.prop_get_wavefront(wf_accessors)),
            "get_amplitude": np.asarray(proper.prop_get_amplitude(wf_accessors)).tolist(),
            "get_phase": np.asarray(proper.prop_get_phase(wf_accessors)).tolist(),
            "after_scalar_add": complex_array_payload(wf_scalar.wfarr),
            "after_matrix_add": complex_array_payload(wf_matrix.wfarr),
        },
        "shape": [n, n],
    }


def main():
    project = Path(__file__).resolve().parents[2]
    pyproper_root = Path(os.environ.get("PYPROPER_ROOT", project / ".." / "proper_v3.3.4_python")).resolve()
    if not (pyproper_root / "proper").is_dir():
        raise RuntimeError(
            f"Missing Python PROPER baseline source tree: {pyproper_root}. "
            "Set PYPROPER_ROOT or place proper_v3.3.4_python next to this repository."
        )
    sys.path.insert(0, str(pyproper_root))

    import proper  # noqa: WPS433

    outdir = Path(__file__).resolve().parent / "baseline" / "python334"
    outdir.mkdir(parents=True, exist_ok=True)

    psf, sampling = run_simple_case(proper)

    csv_path = outdir / "simple_case_psf.csv"
    with csv_path.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerows(psf.tolist())

    meta = {
        "generator": "generate_python_baseline.py",
        "python": sys.version,
        "proper_version": getattr(proper, "__version__", "unknown"),
        "case": "simple_case",
        "shape": list(psf.shape),
        "format": "csv",
        "sampling": sampling,
        "python_root": str(pyproper_root),
        "venv": os.environ.get("VIRTUAL_ENV", ""),
    }
    (outdir / "simple_case_meta.json").write_text(json.dumps(meta, indent=2))

    accessor_case = run_wavefront_accessors_even_case(proper)
    accessor_case.update(
        {
            "generator": "generate_python_baseline.py",
            "python": sys.version,
            "proper_version": getattr(proper, "__version__", "unknown"),
            "case": "wavefront_accessors_even",
            "format": "json",
            "python_root": str(pyproper_root),
            "venv": os.environ.get("VIRTUAL_ENV", ""),
        }
    )
    (outdir / "wavefront_accessors_even.json").write_text(
        json.dumps(accessor_case, indent=2)
    )


if __name__ == "__main__":
    main()
