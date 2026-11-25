"""Download Mesa's lavapipe (Vulkan on the CPU) for Windows: a TEST-ONLY Vulkan driver for CI.

    python tools/get_lavapipe.py --out build/lavapipe [--cache DIR]

Then, for the Vulkan tests (tests/render/test_vulkan_offscreen.cpp):

    $env:PATH = "<OUT>;$env:PATH"                             # vulkan-1.dll (if Windows has none)
    $env:VK_DRIVER_FILES = "<OUT>\\lvp_icd.x86_64.json"         # only this driver
    ctest --test-dir <build> -C Release -L vulkan

In an elevated (administrator) process, which is what CI runners give, the Vulkan loader ignores
VK_DRIVER_FILES and VK_ICD_FILENAMES and reads drivers only from the registry: there, add the
manifest's full path as a DWORD value 0 under HKLM\\SOFTWARE\\Khronos\\Vulkan\\Drivers instead
(.github/workflows/ci.yml does this for its Vulkan test step and removes it afterwards).

Owner decision D14 (docs/specs/m4-plan.md): Windows CI runners have no GPU and no Vulkan driver;
lavapipe stands in for one there, as mesa-vulkan-drivers does in the Linux job. It is never shipped
with MDBridge (users have their graphics driver's Vulkan, or the software renderer).

Source: the MSYS2 project's UCRT64 packages of Mesa (MIT) and of what its vulkan_lvp.dll needs
(LLVM, SPIRV-Tools, the GCC runtime, ...), plus the Khronos Vulkan loader (Apache-2.0), pinned
below by file name and SHA-256 (the hashes of MSYS2's package database, ucrt64.files, on
2026-11-25). Every download is checked against its hash; nothing unchecked is unpacked. OUT gets the
DLLs of NEEDED, the driver manifest lvp_icd.x86_64.json and the packages' licence folders
(licenses/<package>/). The script finally checks that every DLL in OUT finds the DLLs it imports
(in OUT or Windows itself), so a changed dependency fails here, not as a silent test skip.

Updating: take the new package files and SHA256SUM values from MSYS2's database
(https://repo.msys2.org/mingw/ucrt64/ucrt64.files, a tar of <package>/desc files) for Mesa and its
DLL imports, then run this script; its import check names any DLL that is missing.

Needs Python 3.9+ and, before Python 3.14 (compression.zstd), the zstandard module
(pip install zstandard==0.23.0) to unpack .tar.zst.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import struct
import sys
import tarfile
import urllib.request

MIRRORS = ["https://repo.msys2.org/mingw/ucrt64/", "https://mirror.msys2.org/mingw/ucrt64/"]

# (package file, SHA-256): Mesa 26.3.1 and the closure of vulkan_lvp.dll's and vulkan-1.dll's imports
PACKAGES = [
    ("mingw-w64-ucrt-x86_64-mesa-26.3.1-1-any.pkg.tar.zst",
     "8085a7a6acd2e5109c7fbbb8e00781f79a30a4a723197a3fe97d60f8f15be307"),
    ("mingw-w64-ucrt-x86_64-vulkan-loader-1~1.4.363.0-1-any.pkg.tar.zst",
     "065076fe13d39ab97eea4740596426d5ae741b5044e48e1c9495d04835c92dbe"),
    ("mingw-w64-ucrt-x86_64-llvm-libs-23.1.5-1-any.pkg.tar.zst",
     "1739a965f2c3f8f306dfed2575fc30693d3bd0d17fe6e8d6ad99f71b3636e309"),
    ("mingw-w64-ucrt-x86_64-spirv-tools-3~1.4.363.0-1-any.pkg.tar.zst",
     "d31d67fe918aeee21d1e930fd96360856a60146d23094e4ae8baf071a4049aaa"),
    ("mingw-w64-ucrt-x86_64-libstdc++-16.2.0-7-any.pkg.tar.zst",
     "4dd4b216b0e92ffa160127dd50620372609e98cb219ff1c5dcb9b3fbfe153253"),
    ("mingw-w64-ucrt-x86_64-libgcc-16.2.0-7-any.pkg.tar.zst",
     "bc0ae31c582e0847964b4c8c4b55bb2404a9b7ddce7ebf0b176b6693b76e0c62"),
    ("mingw-w64-ucrt-x86_64-libwinpthread-14.0.0.r598.ge3d18511f-1-any.pkg.tar.zst",
     "c35d6be83d141c8cfbb61f1637dd4b53aac8f5562d766aa26d0741c40bad4e4d"),
    ("mingw-w64-ucrt-x86_64-libffi-3.8.1-1-any.pkg.tar.zst",
     "8eeb770c91b7e51ddbcb34a9a991e5b417066270363c3026882149638a70d03d"),
    ("mingw-w64-ucrt-x86_64-libxml2-2.15.4-3-any.pkg.tar.zst",
     "e774c85139c004ce9b4085474d6b3450526c718d2e8c555a87e47adc11f376b1"),
    ("mingw-w64-ucrt-x86_64-zlib-ng-compat-2.3.3-4-any.pkg.tar.zst",
     "477fe283c99e159205ce06dab87610b3ff890dd806289ca1f78b1fdfc6b8f6f8"),
    ("mingw-w64-ucrt-x86_64-zstd-1.5.7-2-any.pkg.tar.zst",
     "dbdb8427280046a2b41697780aa4c52983b708082b0da4755951dc3bea96ca89"),
    ("mingw-w64-ucrt-x86_64-libiconv-1.19.1-1-any.pkg.tar.zst",
     "8a45dc317a285a036bddd79e01505d67b79e44ff49fc3f6e7177a5ea98d098cf"),
    ("mingw-w64-ucrt-x86_64-libsystre-1.0.2-3-any.pkg.tar.zst",
     "62f564b1395d4862db1599b77f6588c32de5c80936845751cf369b31eee1560b"),
    ("mingw-w64-ucrt-x86_64-libtre-0.9.0-2-any.pkg.tar.zst",
     "d50f22545c50d7c0cd584474bba83f0c7a2cd8aab1aed1ad6c1f9c6452f3c4fc"),
    ("mingw-w64-ucrt-x86_64-gettext-runtime-1.0-1-any.pkg.tar.zst",
     "ba693dda4ac375af76ce481ff3a6e7481286546cc7dc6d56c7021dae34084157"),
]

# the DLLs copied to OUT (from the packages' ucrt64/bin)
NEEDED = [
    "vulkan_lvp.dll", "vulkan-1.dll", "libLLVM-23.dll", "libSPIRV-Tools.dll", "libstdc++-6.dll",
    "libgcc_s_seh-1.dll", "libwinpthread-1.dll", "libffi-8.dll", "libxml2-16.dll", "zlib1.dll", "libzstd.dll",
    "libiconv-2.dll", "libsystre-0.dll", "libtre-5.dll", "libintl-8.dll",
]
ICD = "ucrt64/share/vulkan/icd.d/lvp_icd.x86_64.json"

# DLLs Windows provides (imports with these names are not looked for in OUT)
SYSTEM_DLLS = {"advapi32.dll", "cfgmgr32.dll", "gdi32.dll", "kernel32.dll", "ntdll.dll", "ole32.dll",
               "shell32.dll", "user32.dll", "ws2_32.dll", "bcrypt.dll", "dbghelp.dll", "psapi.dll",
               "shlwapi.dll", "version.dll", "msvcrt.dll", "ucrtbase.dll", "oleaut32.dll", "dxgi.dll", "d3d12.dll"}


def zstd_reader(fh):
    """A file object with the decompressed bytes of the zstd stream `fh`."""
    try:
        from compression import zstd  # Python 3.14+
        return zstd.ZstdFile(fh)
    except ImportError:
        pass
    try:
        import zstandard
    except ImportError:
        sys.exit("unpacking .tar.zst needs Python 3.14 or the zstandard module (pip install zstandard==0.23.0)")
    return zstandard.ZstdDecompressor().stream_reader(fh)


def download(name: str, sha256: str, cache: str) -> str:
    path = os.path.join(cache, name)
    if os.path.exists(path) and hashlib.sha256(open(path, "rb").read()).hexdigest() == sha256:
        return path
    errors = []
    for mirror in MIRRORS:
        url = mirror + urllib.request.quote(name)
        try:
            with urllib.request.urlopen(url, timeout=300) as response:
                data = response.read()
        except OSError as exc:
            errors.append(f"{url}: {exc}")
            continue
        digest = hashlib.sha256(data).hexdigest()
        if digest != sha256:
            errors.append(f"{url}: SHA-256 {digest}, expected {sha256}")
            continue
        with open(path + ".part", "wb") as fh:
            fh.write(data)
        os.replace(path + ".part", path)
        print(f"downloaded {name} ({len(data) // 1024} KiB, SHA-256 ok)")
        return path
    sys.exit(f"could not download {name}:\n  " + "\n  ".join(errors) +
             "\n(MSYS2 removes old package versions eventually: update the pins, see the script's help)")


def pe_imports(data: bytes) -> list[str]:
    """Names of the DLLs a PE32+ file imports (import directory only)."""
    pe = struct.unpack_from("<I", data, 0x3C)[0]
    if data[pe:pe + 4] != b"PE\0\0":
        raise ValueError("not a PE file")
    sections_count = struct.unpack_from("<H", data, pe + 6)[0]
    optional_size = struct.unpack_from("<H", data, pe + 20)[0]
    optional = pe + 24
    if struct.unpack_from("<H", data, optional)[0] != 0x20B:
        raise ValueError("not a 64-bit PE file")
    import_rva = struct.unpack_from("<I", data, optional + 112 + 8)[0]  # data directory 1: imports
    sections = []
    for i in range(sections_count):
        at = optional + optional_size + 40 * i
        size, rva, raw_size, raw = struct.unpack_from("<IIII", data, at + 8)
        sections.append((rva, max(size, raw_size), raw))

    def offset(rva: int) -> int:
        for start, size, raw in sections:
            if start <= rva < start + size:
                return raw + rva - start
        raise ValueError(f"RVA {rva:#x} outside the sections")

    names = []
    at = offset(import_rva) if import_rva else None
    while at is not None:
        name_rva = struct.unpack_from("<I", data, at + 12)[0]
        if name_rva == 0:
            break
        start = offset(name_rva)
        names.append(data[start:data.index(b"\0", start)].decode("ascii"))
        at += 20
    return names


def check_imports(out: str) -> None:
    present = {n.lower(): n for n in os.listdir(out) if n.lower().endswith(".dll")}
    missing = []
    for name in sorted(present.values()):
        for dll in pe_imports(open(os.path.join(out, name), "rb").read()):
            low = dll.lower()
            if low in present or low in SYSTEM_DLLS or low.startswith(("api-ms-win-", "ext-ms-win-")):
                continue
            missing.append(f"{name} imports {dll}")
    if missing:
        sys.exit("lavapipe is incomplete (add the packages that provide these DLLs):\n  " + "\n  ".join(missing))


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--out", required=True, help="folder for the DLLs and the driver manifest")
    parser.add_argument("--cache", help="where the downloaded packages are kept (default: OUT/packages)")
    args = parser.parse_args(argv)
    out = os.path.abspath(args.out)
    cache = os.path.abspath(args.cache or os.path.join(out, "packages"))
    os.makedirs(cache, exist_ok=True)

    wanted = {f"ucrt64/bin/{n}": n for n in NEEDED}
    found = {}
    api_version = None
    staging = os.path.join(cache, "bin")  # one package's DLLs at a time
    for name, sha256 in PACKAGES:
        path = download(name, sha256, cache)
        os.makedirs(staging, exist_ok=True)
        links = {}  # hard links within the package (zlib1.dll -> libz-1.dll)
        with open(path, "rb") as fh, tarfile.open(fileobj=zstd_reader(fh), mode="r|") as tar:
            for member in tar:
                if member.islnk() and member.name in wanted:
                    links[member.name] = member.linkname
                    continue
                if not member.isfile():
                    continue
                target = None
                if member.name.startswith("ucrt64/bin/") and member.name.lower().endswith(".dll"):
                    target = os.path.join(staging, os.path.basename(member.name))
                elif member.name.startswith("ucrt64/share/licenses/"):
                    rel = member.name[len("ucrt64/share/licenses/"):]
                    if ".." in rel.split("/"):
                        continue
                    target = os.path.join(out, "licenses", *rel.split("/"))
                elif member.name == ICD:
                    api_version = json.load(tar.extractfile(member))["ICD"]["api_version"]
                    continue
                if target:
                    os.makedirs(os.path.dirname(target), exist_ok=True)
                    with tar.extractfile(member) as src, open(target, "wb") as dst:
                        shutil.copyfileobj(src, dst)
        for member, name in wanted.items():
            source = os.path.join(staging, os.path.basename(links.get(member, member)))
            if member not in found and os.path.exists(source) and (member in links or
                                                                   os.path.basename(member) == name):
                shutil.copyfile(source, os.path.join(out, name))
                found[member] = True
        shutil.rmtree(staging)
    missing = sorted(set(wanted) - set(found))
    if missing or not api_version:
        sys.exit(f"not in the packages: {', '.join(missing) or ICD}")
    # the driver manifest, pointing at vulkan_lvp.dll next to it
    manifest = {"file_format_version": "1.0.1",
                "ICD": {"library_path": ".\\vulkan_lvp.dll", "library_arch": "64", "api_version": api_version}}
    with open(os.path.join(out, "lvp_icd.x86_64.json"), "w", encoding="utf-8", newline="\n") as fh:
        json.dump(manifest, fh, indent=4)
        fh.write("\n")
    check_imports(out)
    print(f"lavapipe (Vulkan {api_version}) in {out}: VK_DRIVER_FILES={os.path.join(out, 'lvp_icd.x86_64.json')}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
