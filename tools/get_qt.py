"""Download an official Qt build: the archives of Qt's online installer, without an account.

    python tools/get_qt.py --version 6.11.3 --arch linux_gcc_64 --out ~/Qt/6.11.3/gcc_64 --archive qtbase --archive icu
    python tools/get_qt.py --version 6.11.3 --arch win64_msvc2022_64 --out C:\\Qt\\6.11.3\\msvc2022_64

OUT becomes the Qt folder (bin, lib, plugins ...) to give CMake as CMAKE_PREFIX_PATH, e.g.
~/Qt/6.11.3/gcc_64. Every archive is checked against the SHA-1 the repository publishes beside it.
Needs py7zr (pip install py7zr).

aqtinstall does the same and more, but its release 3.3.0 cannot read the repository layout Qt uses
for Windows from Qt 6.11 (one folder per compiler); this script reads both layouts.
"""

from __future__ import annotations

import argparse
import hashlib
import os
import sys
import urllib.request
import xml.etree.ElementTree as ET

import py7zr

REPOSITORY = "https://download.qt.io/online/qtsdkrepository"
HOSTS = {"linux": "linux_x64", "win": "windows_x86"}


def fetch(url: str) -> bytes:
    with urllib.request.urlopen(url, timeout=120) as response:
        return response.read()


def find_package(version: str, arch: str):
    """(folder URL, PackageUpdate element) of the Qt build `arch` of `version`."""
    host = HOSTS["win" if arch.startswith("win") else "linux"]
    v = version.replace(".", "")
    base = f"{REPOSITORY}/{host}/desktop/qt6_{v}"
    compiler = arch.split("_", 1)[1]  # win64_msvc2022_64 -> msvc2022_64, linux_gcc_64 -> gcc_64
    name = f"qt.qt6.{v}.{arch}"
    for folder in (f"{base}/qt6_{v}_{compiler}", f"{base}/qt6_{v}"):  # Qt 6.11 Windows layout, then the older one
        try:
            updates = ET.fromstring(fetch(f"{folder}/Updates.xml"))
        except OSError:
            continue
        for package in updates.findall("PackageUpdate"):
            if package.findtext("Name") == name:
                return folder, package
    sys.exit(f"Qt {version} for {arch} not found under {base}")


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--version", required=True, help="e.g. 6.11.3")
    parser.add_argument("--arch", required=True, help="e.g. linux_gcc_64, win64_msvc2022_64")
    parser.add_argument("--out", required=True, help="the Qt folder to create (CMAKE_PREFIX_PATH)")
    parser.add_argument("--archive", action="append", default=[], help="archive to unpack, e.g. qtbase "
                        "(repeatable; default: qtbase)")
    args = parser.parse_args(argv)
    wanted = args.archive or ["qtbase"]
    folder, package = find_package(args.version, args.arch)
    stamp = package.findtext("Version")
    archives = [a.strip() for a in package.findtext("DownloadableArchives", "").split(",")]
    os.makedirs(args.out, exist_ok=True)
    for module in wanted:
        matches = [a for a in archives if a.startswith(module + "-")]
        if not matches:
            sys.exit(f"no '{module}' archive in {package.findtext('Name')}: {', '.join(archives)}")
        for archive in matches:
            url = f"{folder}/{package.findtext('Name')}/{stamp}{archive}"
            data = fetch(url)
            expected = fetch(url + ".sha1").decode().split()[0]
            if hashlib.sha1(data).hexdigest() != expected:
                sys.exit(f"SHA-1 mismatch: {url}")
            target = os.path.join(args.out, archive)
            with open(target, "wb") as fh:
                fh.write(data)
            with py7zr.SevenZipFile(target) as seven:  # ICU's libraries go next to Qt's (Linux)
                seven.extractall(os.path.join(args.out, "lib") if module == "icu" else args.out)
            os.remove(target)
            print(f"{archive} ({len(data) // 1_000_000} MB, sha1 ok)")
    print(f"Qt {args.version} ({args.arch}) in {os.path.abspath(args.out)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
