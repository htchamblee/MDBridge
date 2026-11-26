"""Helpers of the chem generators (make_pdbtools_cases.py, make_annotate_cases.py,
make_builder_cases.py, make_dict_cases.py, make_density_fit_cases.py). Standard library only.

pin_kernels() pins numpy's and OpenBLAS' kernels before numpy is imported, as make_analysis_cases.py
does (OPENBLAS_CORETYPE=Nehalem, one thread, numpy's x86-64 dispatch limited to X86_V2), so that the
helices' 3 x 3 products, Kabsch's SVD, the fits' dot products and the least-squares closures give
the same bits on most machines; regenerate_all.py compares the files whose numbers still depend on
the CPU with the chem tolerance (FITTED there, tests/data/README.md).
"""

import math
import os
import platform


def pin_kernels():
    os.environ["OPENBLAS_CORETYPE"] = "Nehalem"
    for var in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ[var] = "1"
    if platform.machine().lower() in ("x86_64", "amd64"):
        os.environ["NPY_ENABLE_CPU_FEATURES"] = "X86_V2"


def error(exc):
    return {"type": type(exc).__name__, "message": str(exc)}


def fhex(x):
    """float.hex() of a finite float, "nan" / "inf" / "-inf" otherwise."""
    x = float(x)
    if math.isnan(x):
        return "nan"
    if math.isinf(x):
        return "inf" if x > 0 else "-inf"
    return x.hex()


def matrix(m):
    return [[fhex(v) for v in row] for row in m]


def vector(v):
    return [fhex(x) for x in v]
