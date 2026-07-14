#!/usr/bin/env python3
"""Load the pinned Python 3.3.4 baseline with its canonical native kernels."""

import ctypes
import functools
import hashlib
import importlib
import os
import platform
import shlex
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


EXPECTED_PROPER_VERSION = "3.3.4"
KERNEL_SPECS = (
    {
        "name": "cubic_conv",
        "source": "cubic_conv_c.c",
        "library": "libcconv",
        "symbol": "cubic_conv_c",
        "extra_flags": (),
    },
    {
        "name": "cubic_conv_threaded",
        "source": "cubic_conv_threaded_c.c",
        "library": "libcconvthread",
        "symbol": "cubic_conv_c",
        "extra_flags": ("-pthread",),
    },
    {
        "name": "szoom",
        "source": "prop_szoom_c.c",
        "library": "libszoom",
        "symbol": "prop_szoom_c",
        "extra_flags": (),
    },
)


def _sha256_file(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _compiler_command():
    configured = shlex.split(os.environ.get("CC", "cc"))
    if not configured:
        raise RuntimeError("CC is empty; a C compiler is required for Python parity")
    executable = shutil.which(configured[0])
    if executable is None:
        raise RuntimeError(
            f"C compiler {configured[0]!r} was not found; Python parity must use "
            "PROPER's native interpolation kernels"
        )
    return [str(Path(executable).resolve()), *configured[1:]]


def _compiler_version(command):
    completed = subprocess.run(
        [*command, "--version"],
        check=True,
        capture_output=True,
        text=True,
    )
    output = completed.stdout.strip() or completed.stderr.strip()
    if not output:
        raise RuntimeError("C compiler did not report a version")
    return output.splitlines()[0]


class Python334Runtime:
    """Own the temporary shared libraries for one Python baseline process."""

    def __init__(self, pyproper_root):
        self.pyproper_root = Path(pyproper_root).resolve()
        package_root = self.pyproper_root / "proper"
        if not package_root.is_dir():
            raise RuntimeError(
                f"Missing Python PROPER baseline source tree: {self.pyproper_root}"
            )
        if platform.system() not in ("Linux", "Darwin"):
            raise RuntimeError(
                "Pinned Python parity native kernels currently require Linux or macOS"
            )

        self._temporary_directory = tempfile.TemporaryDirectory(
            prefix="proper-python334-kernels-"
        )
        self.build_root = Path(self._temporary_directory.name).resolve()
        compiler = _compiler_command()
        compiler_version = _compiler_version(compiler)
        shared_flag = "-dynamiclib" if platform.system() == "Darwin" else "-shared"
        common_flags = (shared_flag, "-fPIC", "-O2", "-ffp-contract=off", "-fno-fast-math")
        extension = ".dylib" if platform.system() == "Darwin" else ".so"

        libraries = {}
        kernel_metadata = {}
        for spec in KERNEL_SPECS:
            source = package_root / spec["source"]
            if not source.is_file():
                raise RuntimeError(f"Missing required PROPER C source: {source}")
            output = self.build_root / f'{spec["library"]}{extension}'
            flags = (*common_flags, *spec["extra_flags"])
            command = [*compiler, *flags, str(source), "-o", str(output), "-lm"]
            subprocess.run(command, check=True, capture_output=True, text=True)

            library = ctypes.CDLL(str(output))
            if not hasattr(library, spec["symbol"]):
                raise RuntimeError(
                    f'{output} does not export required symbol {spec["symbol"]}'
                )
            libraries[spec["name"]] = output
            kernel_metadata[spec["name"]] = {
                "source": f'proper/{spec["source"]}',
                "source_sha256": _sha256_file(source),
                "compile_flags": [*flags, "-lm"],
                "required_symbol": spec["symbol"],
            }

        if "proper" in sys.modules:
            loaded = Path(sys.modules["proper"].__file__).resolve()
            if self.pyproper_root not in loaded.parents:
                raise RuntimeError(
                    f"A different proper package is already loaded from {loaded}"
                )
            proper = sys.modules["proper"]
        else:
            sys.path.insert(0, str(self.pyproper_root))
            proper = importlib.import_module("proper")

        version = str(getattr(proper, "__version__", "unknown"))
        if version != EXPECTED_PROPER_VERSION:
            raise RuntimeError(
                f"Expected Python PROPER {EXPECTED_PROPER_VERSION}, loaded {version} "
                f"from {Path(proper.__file__).resolve()}"
            )

        proper.cubic_conv_lib = str(libraries["cubic_conv"])
        proper.cubic_conv_threaded_lib = str(libraries["cubic_conv_threaded"])
        proper.szoom_c_lib = str(libraries["szoom"])
        proper.use_cubic_conv = True
        proper.use_cubic_conv_threaded = True
        proper.use_szoom_c = True

        original_cubic_conv = getattr(
            proper.prop_cubic_conv,
            "_proper_native_original",
            proper.prop_cubic_conv,
        )

        @functools.wraps(original_cubic_conv)
        def parity_cubic_conv(
            image_in,
            xval,
            yval,
            THREADED=True,
            GRID=True,
        ):
            coordinate_shape = getattr(xval, "shape", ())
            if not GRID and len(coordinate_shape) == 2:
                if getattr(yval, "shape", ()) != coordinate_shape:
                    raise ValueError("Pointwise X/Y coordinate grids must have equal shapes")
                # The threaded C routine precomputes one X kernel per output
                # column. That is valid for tensor-product grids and separable
                # single-axis tilts, but not when X varies by output row in an
                # arbitrary DM coordinate grid. Use the scalar-equivalent
                # native C path while preserving the upstream X/Y arguments.
                return original_cubic_conv(
                    image_in,
                    xval,
                    yval,
                    THREADED=False,
                    GRID=False,
                )
            return original_cubic_conv(
                image_in,
                xval,
                yval,
                THREADED=THREADED,
                GRID=GRID,
            )

        parity_cubic_conv._proper_native_original = original_cubic_conv
        parity_cubic_conv._proper_dm_coordinate_grid_correction = True
        proper.prop_cubic_conv = parity_cubic_conv
        self.proper = proper
        self.provenance = {
            "required": True,
            "implementation": "upstream PROPER 3.3.4 C sources",
            "platform": platform.system(),
            "machine": platform.machine(),
            "coordinate_grid_policy": (
                "2-D GRID=False uses the serial native C kernel because the "
                "threaded kernel assumes X is constant within each output column"
            ),
            "compiler": {
                "command": compiler,
                "version": compiler_version,
            },
            "kernels": kernel_metadata,
        }
        self.assert_active()

    def assert_active(self):
        expected = {
            "cubic_conv_lib": "cubic_conv",
            "cubic_conv_threaded_lib": "cubic_conv_threaded",
            "szoom_c_lib": "szoom",
        }
        for attribute, kernel_name in expected.items():
            path = Path(getattr(self.proper, attribute)).resolve()
            if path.parent != self.build_root or not path.is_file():
                raise RuntimeError(
                    f"Python parity kernel {kernel_name} is not the pinned temporary build"
                )
        if not all(
            (
                self.proper.use_cubic_conv,
                self.proper.use_cubic_conv_threaded,
                self.proper.use_szoom_c,
            )
        ):
            raise RuntimeError("Python parity attempted to disable a required native kernel")
        if not getattr(
            self.proper.prop_cubic_conv,
            "_proper_dm_coordinate_grid_correction",
            False,
        ):
            raise RuntimeError("Python parity lost its corrected DM coordinate-grid routing")


def load_python334_runtime(pyproper_root):
    """Return a loaded PROPER module and an owner that keeps its kernels alive."""
    runtime = Python334Runtime(pyproper_root)
    return runtime.proper, runtime
