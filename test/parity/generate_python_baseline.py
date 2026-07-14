#!/usr/bin/env python3
import argparse
import csv
import json
import os
import sys
from pathlib import Path

import numpy as np

from baseline_provenance import write_provenance
from python334_runtime import load_python334_runtime


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


def complex_scalar_payload(value):
    value = complex(value)
    return {"real": float(value.real), "imag": float(value.imag)}


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


def run_carrier_phase_case(proper):
    """Capture the upstream opt-in carrier phase used by interferometer arms."""
    n = 4
    wavelength = 500e-9
    prior_phase_offset = getattr(proper, "phase_offset", False)

    def propagated(distance, enabled):
        proper.phase_offset = enabled
        wf = proper.prop_begin(1.0, wavelength, n, beam_diam_fraction=0.5)
        proper.prop_ptp(wf, distance)
        return np.asarray(wf.wfarr)

    try:
        quarter_disabled = propagated(wavelength / 4, False)
        quarter_enabled = propagated(wavelength / 4, True)
        arm1 = propagated(wavelength, True)
        arm2 = propagated(3 * wavelength / 2, True)
    finally:
        proper.phase_offset = prior_phase_offset

    arm_sum = arm1 + arm2
    return {
        "gridsize": n,
        "wavelength_m": wavelength,
        "outputs": {
            "quarter_disabled_mean": complex_scalar_payload(np.mean(quarter_disabled)),
            "quarter_enabled_mean": complex_scalar_payload(np.mean(quarter_enabled)),
            "arm_sum_mean": complex_scalar_payload(np.mean(arm_sum)),
            "arm_mean_intensity": float(np.mean(np.abs(arm_sum) ** 2)),
        },
    }


def run_dm_tilt_projection_case(proper):
    """Exercise asymmetric combined-axis DM projection and public orientation."""
    n = 64
    nact = 8
    yy, xx = np.indices((nact, nact), dtype=np.float64)
    dm_surface = 1e-9 * (
        0.25
        + 0.07 * xx
        - 0.11 * yy
        + 0.013 * xx * yy
        + 0.021 * np.sin((2 * xx + yy) * np.pi / nact)
    )
    config = {
        "beam_diameter_m": 1.0,
        "beam_diam_fraction": 0.5,
        "wavelength_m": 550e-9,
        "grid_size": n,
        "dm_xc": 3.25,
        "dm_yc": 3.75,
        "spacing_m": 0.04,
        "xtilt_deg": 2.3,
        "ytilt_deg": -5.7,
        "ztilt_deg": 1.1,
    }

    def wavefront():
        return proper.prop_begin(
            config["beam_diameter_m"],
            config["wavelength_m"],
            n,
            beam_diam_fraction=config["beam_diam_fraction"],
        )

    kwargs = {
        "XTILT": config["xtilt_deg"],
        "YTILT": config["ytilt_deg"],
        "ZTILT": config["ztilt_deg"],
    }
    wf_zero_tilt = wavefront()
    dmap_zero_tilt = proper.prop_dm(
        wf_zero_tilt,
        dm_surface,
        config["dm_xc"],
        config["dm_yc"],
        config["spacing_m"],
        NO_APPLY=True,
    )
    wf_map = wavefront()
    dmap = proper.prop_dm(
        wf_map,
        dm_surface,
        config["dm_xc"],
        config["dm_yc"],
        config["spacing_m"],
        NO_APPLY=True,
        **kwargs,
    )
    wf_apply = wavefront()
    proper.prop_dm(
        wf_apply,
        dm_surface,
        config["dm_xc"],
        config["dm_yc"],
        config["spacing_m"],
        **kwargs,
    )
    return {
        "config": config,
        "inputs": {"dm_surface_m": np.asarray(dm_surface).tolist()},
        "outputs": {
            "dmap_zero_tilt_m": np.asarray(dmap_zero_tilt).tolist(),
            "dmap_m": np.asarray(dmap).tolist(),
            "centered_wavefront": complex_array_payload(
                proper.prop_get_wavefront(wf_apply)
            ),
        },
        "shape": [n, n],
    }


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
    proper, native_runtime = load_python334_runtime(pyproper_root)
    native_runtime.assert_active()

    outdir = args.output_dir or (
        Path(__file__).resolve().parent / "baseline" / "python334"
    )
    outdir = outdir.resolve()
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

    carrier_case = run_carrier_phase_case(proper)
    carrier_case.update(
        {
            "generator": "generate_python_baseline.py",
            "python": sys.version,
            "proper_version": getattr(proper, "__version__", "unknown"),
            "case": "carrier_phase",
            "format": "json",
            "python_root": str(pyproper_root),
            "venv": os.environ.get("VIRTUAL_ENV", ""),
        }
    )
    (outdir / "carrier_phase.json").write_text(json.dumps(carrier_case, indent=2))

    dm_tilt_case = run_dm_tilt_projection_case(proper)
    dm_tilt_case.update(
        {
            "generator": "generate_python_baseline.py",
            "python": sys.version,
            "proper_version": getattr(proper, "__version__", "unknown"),
            "case": "dm_tilt_projection",
            "format": "json",
            "python_root": str(pyproper_root),
            "venv": os.environ.get("VIRTUAL_ENV", ""),
        }
    )
    (outdir / "dm_tilt_projection.json").write_text(
        json.dumps(dm_tilt_case, indent=2)
    )

    case_metadata = {
        "simple_case": {
            "script": "run_simple_case",
            "wavelength_m": 0.55e-6,
            "grid_size": 256,
            "config": {"beam_diam_fraction": 0.5},
        },
        "wavefront_accessors_even": {
            "script": "run_wavefront_accessors_even_case",
            "wavelength_m": 500e-9,
            "grid_size": 8,
            "config": {"asymmetric_dyadic_inputs": True},
        },
        "carrier_phase": {
            "script": "run_carrier_phase_case",
            "wavelength_m": 500e-9,
            "grid_size": 4,
            "config": {"carrier_phase": [False, True]},
        },
        "dm_tilt_projection": {
            "script": "run_dm_tilt_projection_case",
            "wavelength_m": 550e-9,
            "grid_size": 64,
            "config": {
                "dm_shape": [8, 8],
                "zero_tilt": True,
                "combined_axis_tilts_deg": [2.3, -5.7, 1.1],
                "asymmetric_surface": True,
            },
        },
    }
    write_provenance(
        outdir / "core_baseline_metadata.json",
        generator=Path(__file__).name,
        pyproper_root=pyproper_root,
        proper=proper,
        seed=None,
        cases=case_metadata,
        artifacts=(
            "simple_case_meta.json",
            "simple_case_psf.csv",
            "wavefront_accessors_even.json",
            "carrier_phase.json",
            "dm_tilt_projection.json",
        ),
        native_runtime=native_runtime,
    )


if __name__ == "__main__":
    main()
