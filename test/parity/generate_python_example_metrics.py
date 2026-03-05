#!/usr/bin/env python3
import json
import os
import sys
from pathlib import Path

import numpy as np


def summarize(a, sampling):
    arr = np.asarray(a)
    out = {
        "shape": list(arr.shape),
        "sampling": float(sampling),
    }
    if np.iscomplexobj(arr):
        out.update(
            {
                "sum_re": float(np.sum(arr.real)),
                "sum_im": float(np.sum(arr.imag)),
                "norm": float(np.linalg.norm(arr)),
                "center_re": float(arr[arr.shape[0] // 2, arr.shape[1] // 2].real),
                "center_im": float(arr[arr.shape[0] // 2, arr.shape[1] // 2].imag),
            }
        )
    else:
        out.update(
            {
                "sum": float(np.sum(arr)),
                "norm": float(np.linalg.norm(arr)),
                "max": float(np.max(arr)),
                "center": float(arr[arr.shape[0] // 2, arr.shape[1] // 2]),
            }
        )
    return out


def main():
    project = Path(__file__).resolve().parents[2]
    pyproper_root = (project / ".." / "proper_v3.3.4_python").resolve()
    pyexamples = pyproper_root / "proper" / "examples"
    sys.path.insert(0, str(pyproper_root))
    sys.path.insert(0, str(pyexamples))

    import proper
    import simple_prescription
    import simple_telescope
    import hubble_simple
    import microscope
    import example_system
    import psdtest
    import talbot
    import talbot_correct
    import run_occulter
    import run_coronagraph
    import run_coronagraph_dm
    import multi_example

    np.random.seed(12345)

    cases = {
        "simple_prescription": simple_prescription.simple_prescription(0.55e-6, 256),
        "simple_telescope": simple_telescope.simple_telescope(0.55e-6, 256),
        "hubble_simple": hubble_simple.hubble_simple(0.55e-6, 256, PASSVALUE={"delta_sec": 0.0}),
        "microscope": microscope.microscope(0.55e-6, 256, PASSVALUE={"focus_offset": 0.0}),
        "example_system": example_system.example_system(0.55e-6, 256),
        "psdtest": psdtest.psdtest(0.55e-6, 256, PASSVALUE={"usepsdmap": True}),
        "talbot": talbot.talbot(0.5e-6, 128, PASSVALUE={"diam": 0.1, "period": 0.04, "dist": 0.0}),
        "talbot_correct": talbot_correct.talbot_correct(0.5e-6, 128, PASSVALUE={"diam": 0.1, "period": 0.04, "dist": 0.0}),
        "run_occulter": run_occulter.run_occulter(0.55e-6, 256, PASSVALUE={"occulter_type": "GAUSSIAN"}),
        "run_coronagraph": run_coronagraph.run_coronagraph(0.55e-6, 256, PASSVALUE={"use_errors": False, "occulter_type": "GAUSSIAN"}),
        "run_coronagraph_dm": run_coronagraph_dm.run_coronagraph_dm(0.55e-6, 256, PASSVALUE={"use_errors": False, "use_dm": False, "occulter_type": "GAUSSIAN"}),
        "multi_example": multi_example.multi_example(0.55e-6, 256, PASSVALUE={"use_dm": False, "dm": np.zeros((48, 48), dtype=np.float64)}),
    }

    metrics = {name: summarize(arr, sampling) for name, (arr, sampling) in cases.items()}

    outdir = Path(__file__).resolve().parent / "baseline" / "python334"
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / "example_metrics.json").write_text(json.dumps(metrics, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
