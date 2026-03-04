#!/usr/bin/env python3
import csv
import json
from pathlib import Path


def main():
    outdir = Path(__file__).resolve().parent / "baseline" / "python334"
    outdir.mkdir(parents=True, exist_ok=True)

    n = 256
    csv_path = outdir / "simple_case_psf.csv"
    with csv_path.open("w", newline="") as f:
        w = csv.writer(f)
        row = [1.0] * n
        for _ in range(n):
            w.writerow(row)

    meta = {
        "generator": "generate_python_baseline.py",
        "case": "simple_case",
        "shape": [n, n],
        "format": "csv",
        "note": "placeholder baseline",
    }
    (outdir / "simple_case_meta.json").write_text(json.dumps(meta, indent=2))


if __name__ == "__main__":
    main()
