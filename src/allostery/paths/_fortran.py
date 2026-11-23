"""Loading (and, if needed, compiling) the Fortran path-sampling kernel.

The kernel ``_twisp.f90`` is a small, dependency-free Fortran file that is
compiled into a shared library and called through :mod:`ctypes`.  The library
is looked up in this order:

1. the file named by the environment variable ``ALLOSTERY_TWISP_LIB``;
2. a library next to this module (built by ``allostery compile``);
3. a per-user cache directory (``$ALLOSTERY_CACHE_DIR`` or
   ``~/.cache/allostery``), where the library is compiled automatically the
   first time it is needed, using ``$FC`` or ``gfortran``.

If none of this works, :func:`load` returns ``None`` and the (much slower)
pure-Python implementation is used instead.
"""

from __future__ import annotations

import ctypes
import hashlib
import os
import shutil
import subprocess
import sys
import tempfile
import threading
import warnings
from pathlib import Path

import numpy as np

SOURCE = Path(__file__).with_name("_twisp.f90")
_LIB_BASENAME = "libtwisp"
if sys.platform == "darwin":
    _EXT = ".dylib"
elif sys.platform.startswith("win"):
    _EXT = ".dll"
else:
    _EXT = ".so"

DEFAULT_FLAGS = ["-O2", "-fPIC", "-shared"]

_lock = threading.Lock()
_lib = None
_load_attempted = False
_load_error: str | None = None


def _source_hash() -> str:
    return hashlib.sha256(SOURCE.read_bytes()).hexdigest()[:12]


def _cache_dir() -> Path:
    root = os.environ.get("ALLOSTERY_CACHE_DIR")
    if root:
        return Path(root)
    xdg = os.environ.get("XDG_CACHE_HOME")
    return Path(xdg) / "allostery" if xdg else Path.home() / ".cache" / "allostery"


def find_compiler() -> str | None:
    """Return the Fortran compiler to use, or ``None``."""
    fc = os.environ.get("FC")
    if fc and shutil.which(fc.split()[0]):
        return fc
    for name in ("gfortran", "flang-new", "flang", "ifx", "ifort"):
        if shutil.which(name):
            return name
    return None


def compile_library(output: str | os.PathLike | None = None, compiler: str | None = None,
                    flags: list[str] | None = None, verbose: bool = False) -> Path:
    """Compile the Fortran kernel into a shared library and return its path.

    Parameters
    ----------
    output
        Destination file.  Defaults to the per-user cache location.
    compiler
        Fortran compiler (default: ``$FC`` or ``gfortran``).
    flags
        Compiler flags (default ``-O2 -fPIC -shared``).  Do not use
        ``-ffast-math``: the algorithm relies on IEEE Inf/NaN semantics.
    """
    compiler = compiler or find_compiler()
    if compiler is None:
        raise RuntimeError("No Fortran compiler found (install gfortran or set $FC).")
    if output is None:
        output = _cache_dir() / library_name()
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    if flags is None:
        flags = list(DEFAULT_FLAGS)
        if "gfortran" in compiler:
            # no fused multiply-add: keeps results identical to the
            # pure-Python implementation on every platform
            flags.append("-ffp-contract=off")
    flags = list(flags)
    if sys.platform == "darwin" and "-shared" in flags:
        flags = [f for f in flags if f != "-shared"] + ["-dynamiclib"]
    # compile into a private temporary directory and atomically move the
    # result into place, so that concurrent processes do not collide.
    with tempfile.TemporaryDirectory() as tmp:
        tmp_lib = Path(tmp) / output.name
        cmd = compiler.split() + flags + ["-J", tmp, str(SOURCE), "-o", str(tmp_lib)]
        if verbose:
            print(" ".join(cmd))
        res = subprocess.run(cmd, capture_output=True, text=True)
        if res.returncode != 0:
            raise RuntimeError(f"Compilation of {SOURCE.name} failed:\n{res.stderr}")
        tmp_dest = output.with_name(output.name + f".tmp{os.getpid()}")
        shutil.copy2(tmp_lib, tmp_dest)
        os.replace(tmp_dest, output)
    return output


def _candidates():
    env = os.environ.get("ALLOSTERY_TWISP_LIB")
    if env:
        yield Path(env)
    yield package_library_path()
    yield _cache_dir() / library_name()


def library_name() -> str:
    """File name of the library for the current kernel source.

    The name contains a hash of ``_twisp.f90``, so a library compiled from an
    older version of the source is never loaded by mistake.
    """
    return f"{_LIB_BASENAME}-{_source_hash()}{_EXT}"


def package_library_path() -> Path:
    return Path(__file__).with_name(library_name())


def _configure(lib):
    dptr = np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS")
    iptr = np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS")
    c_int, c_double = ctypes.c_int, ctypes.c_double
    pint = ctypes.POINTER(c_int)
    pdouble = ctypes.POINTER(c_double)

    lib.allostery_twisp_dijkstra.restype = None
    lib.allostery_twisp_dijkstra.argtypes = [c_int, dptr, c_int, c_int, iptr, pint, pdouble]

    lib.allostery_twisp_sample.restype = None
    lib.allostery_twisp_sample.argtypes = [
        c_int, dptr, dptr, c_int, c_int, c_int, c_double, c_int,
        iptr, pint, iptr, c_int, iptr, iptr, dptr, pint, pint, pint,
    ]
    return lib


def load(compile_if_missing: bool = True):
    """Return the loaded kernel (a :class:`ctypes.CDLL`) or ``None``."""
    global _lib, _load_attempted, _load_error
    with _lock:
        if _lib is not None or _load_attempted:
            return _lib
        _load_attempted = True
        if os.environ.get("ALLOSTERY_DISABLE_FORTRAN"):
            _load_error = "disabled through ALLOSTERY_DISABLE_FORTRAN"
            return None
        for path in _candidates():
            if path.is_file():
                try:
                    _lib = _configure(ctypes.CDLL(str(path)))
                    return _lib
                except OSError as exc:  # pragma: no cover - broken library
                    _load_error = str(exc)
        if compile_if_missing:
            try:
                path = compile_library()
                _lib = _configure(ctypes.CDLL(str(path)))
                return _lib
            except Exception as exc:  # pragma: no cover - no compiler
                _load_error = str(exc)
        warnings.warn(
            "The Fortran path-sampling kernel is not available "
            f"({_load_error or 'no compiled library found'}); falling back to the "
            "slow pure-Python implementation.  Install gfortran and run "
            "`allostery compile` to fix this.",
            RuntimeWarning,
            stacklevel=3,
        )
        return None


def available() -> bool:
    return load() is not None


def dijkstra(cost: np.ndarray, source: int, sink: int):
    """Fortran Dijkstra on the edge costs of :func:`._reference.cost_matrix`;
    returns (path as 1-based array, length)."""
    lib = load()
    n = cost.shape[0]
    path = np.zeros(n, dtype=np.int32)
    npath = ctypes.c_int(0)
    dist = ctypes.c_double(0.0)
    lib.allostery_twisp_dijkstra(n, cost, source + 1, sink + 1, path,
                                 ctypes.byref(npath), ctypes.byref(dist))
    return path[: npath.value].copy(), dist.value


def sample(cost, positions, angular, source, sink, temperature, nsweeps,
           path, npath, state, nodes_buf, lengths_buf, weights_buf):
    """Run the Fortran sampler for at most ``nsweeps`` paths.

    ``path``/``state`` are modified in place.  Returns
    ``(new_npath, ndone, nnodes, status)``.
    """
    lib = load()
    n = cost.shape[0]
    c_npath = ctypes.c_int(npath)
    ndone = ctypes.c_int(0)
    nnodes = ctypes.c_int(0)
    status = ctypes.c_int(0)
    lib.allostery_twisp_sample(
        n, cost, positions, int(angular), source + 1, sink + 1,
        float(temperature), int(nsweeps), path, ctypes.byref(c_npath), state,
        nodes_buf.size, nodes_buf, lengths_buf, weights_buf,
        ctypes.byref(ndone), ctypes.byref(nnodes), ctypes.byref(status),
    )
    return c_npath.value, ndone.value, nnodes.value, status.value
