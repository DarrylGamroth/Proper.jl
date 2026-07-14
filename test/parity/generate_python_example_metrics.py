#!/usr/bin/env python3
import argparse
import json
import os
import sys
from pathlib import Path

import numpy as np

from baseline_provenance import write_provenance
from python334_runtime import load_python334_runtime


def probe_indices(shape):
    ny, nx = shape
    return (
        (ny // 5, nx // 3),
        (ny // 3, 2 * nx // 3),
        (3 * ny // 5, nx // 4),
        (4 * ny // 5, 3 * nx // 4),
    )


def summarize(a, sampling):
    arr = np.asarray(a)
    probes = [arr[index] for index in probe_indices(arr.shape)]
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
                "probes_re": [float(value.real) for value in probes],
                "probes_im": [float(value.imag) for value in probes],
            }
        )
    else:
        out.update(
            {
                "sum": float(np.sum(arr)),
                "norm": float(np.linalg.norm(arr)),
                "max": float(np.max(arr)),
                "center": float(arr[arr.shape[0] // 2, arr.shape[1] // 2]),
                "probes": [float(value) for value in probes],
            }
        )
    return out


def zero_unwritten_szoom_border(image, input_size, magnification):
    """Clear pixels skipped by prop_szoom_c when its complex output used np.empty."""
    output_size = image.shape[0]

    def nearest(value):
        return np.floor(value + 0.5) if value > 0 else -np.floor(-value + 0.5)

    valid = np.zeros(output_size, dtype=bool)
    for output_index in range(output_size):
        input_offset = (output_index - output_size // 2) / magnification
        input_pixel = int(nearest(input_offset)) + input_size // 2
        valid[output_index] = input_pixel - 6 >= 0 and input_pixel + 6 < input_size
    image[~np.outer(valid, valid)] = 0
    return image


def run_testmulti1(proper):
    lambda_min = 0.5
    lambda_max = 0.7
    nlambda = 9
    gridsize = 256
    npsf = 256
    final_sampling = 1.5e-6

    wavelength = (
        np.arange(nlambda) / (nlambda - 1.0) * (lambda_max - lambda_min)
        + lambda_min
    )
    passvalue = {"use_dm": True, "dm": np.zeros((48, 48), dtype=np.float64)}
    passvalue["dm"][20, 20] = 0.2e-6
    passvalue["dm"][15, 25] = 0.2e-6

    fields, sampling = proper.prop_run_multi(
        "multi_example", wavelength, gridsize, PASSVALUE=passvalue
    )
    psfs = np.zeros((nlambda, npsf, npsf), dtype=np.float64)
    for index in range(nlambda):
        magnification = sampling[index] / final_sampling
        field = proper.prop_magnify(
            fields[index, :, :], magnification, npsf, CONSERVE=True
        )
        zero_unwritten_szoom_border(field, fields.shape[-1], magnification)
        psfs[index, :, :] = np.abs(field) ** 2
    return np.mean(psfs, axis=0), final_sampling


def run_testmulti2(proper):
    wavelength = 0.6
    gridsize = 256
    npatterns = 3
    passvalues = [
        {"use_dm": True, "dm": np.zeros((48, 48), dtype=np.float64)}
        for _ in range(npatterns)
    ]
    x = np.dot(
        (np.arange(48.0) / 47 * (2 * np.pi)).reshape(48, 1),
        np.ones((1, 48), dtype=np.float64),
    )
    for index in range(npatterns):
        passvalues[index]["dm"] = 5.0e-8 * np.cos(4 * x * (index + 1))
    return proper.prop_run_multi(
        "multi_example", wavelength, gridsize, PASSVALUE=passvalues
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args()

    project = Path(__file__).resolve().parents[2]
    pyproper_root = Path(os.environ.get("PYPROPER_ROOT", project / ".." / "proper_v3.3.4_python")).resolve()
    if not (pyproper_root / "proper").is_dir():
        raise RuntimeError(
            f"Missing Python PROPER baseline source tree: {pyproper_root}. "
            "Set PYPROPER_ROOT or place proper_v3.3.4_python next to this repository."
        )
    pyexamples = pyproper_root / "proper" / "examples"
    sys.path.insert(0, str(pyexamples))

    proper, native_runtime = load_python334_runtime(pyproper_root)
    native_runtime.assert_active()
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

    cases["testmulti1"] = run_testmulti1(proper)
    testmulti2_fields, testmulti2_sampling = run_testmulti2(proper)
    for index in range(testmulti2_fields.shape[0]):
        cases[f"testmulti2_pattern_{index + 1}"] = (
            testmulti2_fields[index, :, :],
            testmulti2_sampling[index],
        )

    metrics = {name: summarize(arr, sampling) for name, (arr, sampling) in cases.items()}

    outdir = args.output_dir or (
        Path(__file__).resolve().parent / "baseline" / "python334"
    )
    outdir = outdir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    metrics_path = outdir / "example_metrics.json"
    metrics_path.write_text(json.dumps(metrics, indent=2, sort_keys=True) + "\n")

    case_metadata = {
        "simple_prescription": {"wavelength_m": 0.55e-6, "grid_size": 256, "config": {}},
        "simple_telescope": {"wavelength_m": 0.55e-6, "grid_size": 256, "config": {}},
        "hubble_simple": {"wavelength_m": 0.55e-6, "grid_size": 256, "config": {"delta_sec": 0.0}},
        "microscope": {"wavelength_m": 0.55e-6, "grid_size": 256, "config": {"focus_offset": 0.0}},
        "example_system": {"wavelength_m": 0.55e-6, "grid_size": 256, "config": {}},
        "psdtest": {"wavelength_m": 0.55e-6, "grid_size": 256, "config": {"usepsdmap": True}},
        "talbot": {"wavelength_m": 0.5e-6, "grid_size": 128, "config": {"diam": 0.1, "period": 0.04, "dist": 0.0}},
        "talbot_correct": {"wavelength_m": 0.5e-6, "grid_size": 128, "config": {"diam": 0.1, "period": 0.04, "dist": 0.0}},
        "run_occulter": {"wavelength_m": 0.55e-6, "grid_size": 256, "config": {"occulter_type": "GAUSSIAN"}},
        "run_coronagraph": {"wavelength_m": 0.55e-6, "grid_size": 256, "config": {"use_errors": False, "occulter_type": "GAUSSIAN"}},
        "run_coronagraph_dm": {"wavelength_m": 0.55e-6, "grid_size": 256, "config": {"use_errors": False, "use_dm": False, "occulter_type": "GAUSSIAN"}},
        "multi_example": {"wavelength_m": 0.55e-6, "grid_size": 256, "config": {"use_dm": False}},
        "testmulti1": {"wavelength_m": [0.5e-6, 0.7e-6], "grid_size": 256, "config": {"nlambda": 9, "dm_pokes_m": [0.2e-6, 0.2e-6], "final_sampling_m": 1.5e-6, "zero_unwritten_szoom_c_border": True}},
        "testmulti2_pattern_1": {"wavelength_m": 0.6e-6, "grid_size": 256, "config": {"pattern": 1, "amplitude_m": 5e-8}},
        "testmulti2_pattern_2": {"wavelength_m": 0.6e-6, "grid_size": 256, "config": {"pattern": 2, "amplitude_m": 5e-8}},
        "testmulti2_pattern_3": {"wavelength_m": 0.6e-6, "grid_size": 256, "config": {"pattern": 3, "amplitude_m": 5e-8}},
    }
    write_provenance(
        outdir / "example_metrics_metadata.json",
        generator=Path(__file__).name,
        pyproper_root=pyproper_root,
        proper=proper,
        seed=12345,
        cases=case_metadata,
        artifacts=("example_metrics.json",),
        native_runtime=native_runtime,
    )


if __name__ == "__main__":
    main()
