"""Copy the FFT that SciPy 1.18.1 runs (ducc0 0.40.0, as SciPy vendors it) into third_party/ducc0.

    python -I -B tools/vendor_ducc0.py SCIPY_SDIST [--check]

SCIPY_SDIST is SciPy's source release scipy-1.18.1.tar.gz (PyPI), checked against SDIST_SHA256 below.
`scipy.fft` is built on Martin Reinecke's ducc0 FFT (the successor of pocketfft): the sdist holds the
subset SciPy compiles in subprojects/duccfft/ducc0/ (tools/vendoring/vendor_duccfft.sh of SciPy's
repository, ducc0 0.40.0). Owner decision D3 (docs/specs/m6-plan.md) vendors that code so that
MDBridge's rfftn / irfftn give SciPy's numbers (src/numerics/fft.cpp).

The files are copied unmodified, with ducc0's licence (BSD-3-Clause; the sources say
"SPDX-License-Identifier: BSD-3-Clause OR GPL-2.0-or-later", MDBridge takes them under BSD-3-Clause),
except SciPy's pybind11 bindings (ducc0/bindings) and ducc0/math/pswf_rokhlin.h, which the FFT does not
use. third_party/ducc0/SHA256SUMS lists every file's SHA-256; src/numerics/CMakeLists.txt checks them
at configure time, so an edited copy fails the build. --check compares third_party/ducc0 with the
sdist instead of writing it. Standard library only.
"""

import argparse
import hashlib
import os
import sys
import tarfile

SDIST_SHA256 = "52c4b7422442aba924d03ad4019852b08a92e64ea187b933135687bfe2747307"  # scipy-1.18.1.tar.gz
PREFIX = "scipy-1.18.1/subprojects/duccfft/"
FILES = [
    "LICENSE.md",
    "ducc0/fft/fft.h",
    "ducc0/fft/fft1d_impl.h",
    "ducc0/fft/fftnd_impl.h",
    "ducc0/infra/aligned_array.h",
    "ducc0/infra/error_handling.h",
    "ducc0/infra/mav.cc",
    "ducc0/infra/mav.h",
    "ducc0/infra/misc_utils.h",
    "ducc0/infra/simd.h",
    "ducc0/infra/string_utils.cc",
    "ducc0/infra/string_utils.h",
    "ducc0/infra/threading.cc",
    "ducc0/infra/threading.h",
    "ducc0/infra/useful_macros.h",
    "ducc0/math/cmplx.h",
    "ducc0/math/unity_roots.h",
]
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "third_party", "ducc0")


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("sdist")
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    with open(args.sdist, "rb") as f:
        digest = hashlib.sha256(f.read()).hexdigest()
    if digest != SDIST_SHA256:
        sys.exit(f"error: {args.sdist}: SHA-256 {digest}, expected {SDIST_SHA256} (scipy-1.18.1.tar.gz)")
    sums = []
    with tarfile.open(args.sdist, "r:gz") as tar:
        for rel in FILES:
            data = tar.extractfile(PREFIX + rel).read()
            sums.append(f"{hashlib.sha256(data).hexdigest()}  {rel}\n")
            path = os.path.join(OUT, rel)
            if args.check:
                with open(path, "rb") as f:
                    if f.read() != data:
                        sys.exit(f"error: {path} differs from SciPy 1.18.1's copy")
            else:
                os.makedirs(os.path.dirname(path), exist_ok=True)
                with open(path, "wb") as f:
                    f.write(data)
    text = "".join(sums)
    sums_path = os.path.join(OUT, "SHA256SUMS")
    if args.check:
        with open(sums_path, encoding="utf-8") as f:
            if f.read() != text:
                sys.exit(f"error: {sums_path} is out of date")
        print(f"third_party/ducc0: {len(FILES)} files equal to SciPy 1.18.1's")
    else:
        with open(sums_path, "w", encoding="utf-8", newline="\n") as f:
            f.write(text)
        print(f"wrote {len(FILES)} files and SHA256SUMS to {OUT}")


if __name__ == "__main__":
    main()
