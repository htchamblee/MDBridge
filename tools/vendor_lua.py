"""Copy Lua 5.4.9 (PUC-Rio's release) into third_party/lua.

    python -I -B tools/vendor_lua.py LUA_TARBALL [--check]

LUA_TARBALL is the official release lua-5.4.9.tar.gz (https://www.lua.org/ftp/lua-5.4.9.tar.gz,
released 2026-08-10), checked against TARBALL_SHA256 below (the SHA-256 lua.org lists for it).
Milestone 7's owner decision D1 (docs/specs/m7-plan.md) vendors Lua 5.4 for the workflow generators:
the interpreter and the libraries the sandbox opens are compiled as C++ by src/workflow/lua.

The files are copied unmodified: the library sources (src/, without the stand-alone programs lua.c
and luac.c and the Makefile), the README and doc/readme.html, which holds the MIT licence (also at the
end of src/lua.h). third_party/lua/SHA256SUMS lists every file's SHA-256; src/workflow/lua/CMakeLists.txt
checks them at configure time, so an edited copy fails the build. --check compares third_party/lua
with the release instead of writing it. Standard library only.
"""

import argparse
import hashlib
import os
import sys
import tarfile

VERSION = "5.4.9"
TARBALL_SHA256 = "2335b6c582a52654f94612bf10d2f4672805d05329aa6568b1d8cd9e5c6fb8e6"  # lua-5.4.9.tar.gz
PREFIX = f"lua-{VERSION}/"
EXCLUDED = {"src/Makefile", "src/lua.c", "src/luac.c"}  # the stand-alone programs and their build
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "third_party", "lua")


def release_files(tar):
    """The files to vendor, in sorted order: README, doc/readme.html and src/ without EXCLUDED."""
    names = []
    for member in tar.getmembers():
        if not member.isfile() or not member.name.startswith(PREFIX):
            continue
        rel = member.name[len(PREFIX):]
        if rel in ("README", "doc/readme.html") or (rel.startswith("src/") and rel not in EXCLUDED):
            names.append(rel)
    return sorted(names)


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("tarball")
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    with open(args.tarball, "rb") as f:
        digest = hashlib.sha256(f.read()).hexdigest()
    if digest != TARBALL_SHA256:
        sys.exit(f"error: {args.tarball}: SHA-256 {digest}, expected {TARBALL_SHA256} (lua-{VERSION}.tar.gz)")
    sums = []
    with tarfile.open(args.tarball, "r:gz") as tar:
        files = release_files(tar)
        for rel in files:
            data = tar.extractfile(PREFIX + rel).read()
            sums.append(f"{hashlib.sha256(data).hexdigest()}  {rel}\n")
            path = os.path.join(OUT, *rel.split("/"))
            if args.check:
                if not os.path.isfile(path):
                    sys.exit(f"error: {path} is missing")
                with open(path, "rb") as f:
                    if f.read() != data:
                        sys.exit(f"error: {path} differs from Lua {VERSION}'s copy")
            else:
                os.makedirs(os.path.dirname(path), exist_ok=True)
                with open(path, "wb") as f:
                    f.write(data)
    if args.check:  # nothing else may sit beside the release's files
        known = set(files) | {"SHA256SUMS"}
        for folder, _, names in os.walk(OUT):
            for name in names:
                rel = os.path.relpath(os.path.join(folder, name), OUT).replace(os.sep, "/")
                if rel not in known:
                    sys.exit(f"error: third_party/lua/{rel} is not part of Lua {VERSION}'s release")
    text = "".join(sums)
    sums_path = os.path.join(OUT, "SHA256SUMS")
    if args.check:
        with open(sums_path, encoding="utf-8") as f:
            if f.read() != text:
                sys.exit(f"error: {sums_path} is out of date")
        print(f"third_party/lua: {len(files)} files equal to Lua {VERSION}'s")
    else:
        with open(sums_path, "w", encoding="utf-8", newline="\n") as f:
            f.write(text)
        print(f"wrote {len(files)} files and SHA256SUMS to {OUT}")


if __name__ == "__main__":
    main()
