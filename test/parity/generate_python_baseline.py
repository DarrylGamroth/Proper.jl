#!/usr/bin/env python3
import csv
import json
import os
import sys
from pathlib import Path


def run_simple_case(proper):
    wf = proper.prop_begin(2.4, 0.55e-6, 256, beam_diam_fraction=0.5)
    proper.prop_circular_aperture(wf, 0.6)
    proper.prop_lens(wf, 20.0)
    proper.prop_propagate(wf, 20.0)
    psf, sampling = proper.prop_end(wf)
    return psf, float(sampling)


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


if __name__ == "__main__":
    main()
