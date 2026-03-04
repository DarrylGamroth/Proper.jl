#!/usr/bin/env python3
import json
import numpy as np
from pathlib import Path


def main():
    outdir = Path(__file__).resolve().parent / "baseline" / "python334"
    outdir.mkdir(parents=True, exist_ok=True)

    # Placeholder baseline until full prescription parity harness lands.
    psf = np.ones((256, 256), dtype=np.float64)
    np.save(outdir / "simple_case_psf.npy", psf)
    meta = {
        "generator": "generate_python_baseline.py",
        "case": "simple_case",
        "shape": list(psf.shape),
        "note": "placeholder baseline",
    }
    (outdir / "simple_case_meta.json").write_text(json.dumps(meta, indent=2))


if __name__ == "__main__":
    main()
