"""Record what MDATK reads from project files: the reference for MDBridge's reader.

    python tests/parity/reference.py --out DIR FILE...

For each FILE, MDATK's mdatk.projectfile (v0.0.4.1) gives the expected result, written to DIR as

* ``<name>.settings.json``  projectfile.read(FILE)
* ``<name>.log.json``       projectfile.read_log(FILE)
* ``<name>.error.txt``      instead of both, when read() rejects the file (the message, with the
                            file's path replaced by its name)

JSON text is json.dumps(..., indent=2, ensure_ascii=False) with dates as isoformat(), which
MDBridge's core::toJson() reproduces byte for byte.
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import logging
import os
import sys

from mdatk import projectfile as PF
from mdatk import tomlio


def dumps(data) -> str:
    def default(value):
        if isinstance(value, (dt.datetime, dt.date, dt.time)):
            return value.isoformat()
        raise TypeError(f"cannot write {type(value).__name__}")

    return json.dumps(data, indent=2, ensure_ascii=False, default=default) + "\n"


def record(path: str, out: str) -> str:
    name = os.path.basename(path)
    for suffix in (".settings.json", ".log.json", ".error.txt"):
        target = os.path.join(out, name + suffix)
        if os.path.exists(target):
            os.remove(target)
    try:
        settings = PF.read(path)
    except (ValueError, tomlio.TOMLDecodeError) as exc:
        with open(os.path.join(out, name + ".error.txt"), "w", encoding="utf-8") as fh:
            fh.write(str(exc).replace(path, name) + "\n")
        return "error"
    with open(os.path.join(out, name + ".settings.json"), "w", encoding="utf-8") as fh:
        fh.write(dumps(settings))
    with open(os.path.join(out, name + ".log.json"), "w", encoding="utf-8") as fh:
        fh.write(dumps(PF.read_log(path)))
    return "ok"


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--out", required=True, help="folder for the reference files")
    parser.add_argument("files", nargs="+")
    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.ERROR)
    os.makedirs(args.out, exist_ok=True)
    for path in args.files:
        print(f"{record(path, args.out):5s} {path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
