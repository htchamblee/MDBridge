"""What MDATK's map reader and writer and the map functions of mdatk/density.py give
(tests/data/io/maps, for tests/io/test_mrc.cpp).

    /opt/mdbridge-venv/bin/python -I -B make_map_cases.py MDATK_CHECKOUT OUT_DIR

Needs numpy and scipy (mdatk.density imports them); nothing is written into the checkout. Writes

  files/*             small MRC / CCP4 files, plain and gzip: every axis order in both byte
                      orders, modes 0 (signed / unsigned rules), 1, 2, 6, 12, machine stamps and
                      byte-order guessing, NSTART / ORIGIN / MXYZ / cell / angle rules, extended
                      headers, gzip members, flags, padding, truncation and corruption, every
                      MapError
  read_cases.json     mrc.read_map of each file: grid, origin, voxel, SHA-256 of the data (float32
                      little-endian, [z][y][x]), stats() and the first values; or the error. "error"
                      is MDBridge's message (MDATK's for MapError, the documented deviation for
                      Python's EOFError / zlib.error / BadGzipFile / MemoryError, which "python"
                      records)
  info_cases.json     the output of `mdatk density info` for some of the files
  stats_cases.json    DensityMap.stats() (numpy's float32 min / max / mean / std, as float32 bit
                      patterns) and map_level of generated float32 arrays (below), up to 2^24 + 1
                      values (numpy's pairwise summation)
  write_cases.json    write_map of generated maps: SHA-256 of the file, its first 224 header bytes
                      (the rest is zero) and what read_map reads back
  simulate_cases.json atom_weights, and simulate_map of atom sets: grid, origin, voxel, SHA-256 of
                      the data, stats. numpy's exp is replaced by the C library's (math.exp, which
                      numpy uses without SIMD kernels); "numpy_sha256" / "numpy_differing_voxels"
                      record what numpy's own exp gives on the generating machine

Generated arrays ("data" specs), reproduced by the C++ test: z_i = splitmix64(seed + i * golden),
i = 1..n, then
  uniform   float32((((z >> 40) / 2^24) - 0.5) * scale + offset)       (float64 arithmetic)
  bits      float32 bits (z & 0x807FFFFF) | (emin + ((z >> 32) % (emax - emin + 1))) << 23
  constant  float32(value)
  ramp      float32(start + i0 * step), i0 = 0..n-1
Floats are JSON numbers (NaN / Infinity as json writes them); file paths in messages are "{path}".
"""

import contextlib
import gzip
import hashlib
import io
import json
import math
import os
import resource
import shutil
import struct
import sys
import warnings
import zlib

warnings.filterwarnings("ignore")

if len(sys.argv) != 3:
    sys.exit("usage: make_map_cases.py MDATK_CHECKOUT OUT_DIR")
sys.path.insert(0, os.path.abspath(sys.argv[1]))

import numpy as np  # noqa: E402

from mdatk import cli  # noqa: E402
from mdatk import density as D  # noqa: E402
from mdatk.io.mrc import DensityMap, MapError, read_map, write_map  # noqa: E402

OUT = os.path.abspath(sys.argv[2])
FILES = os.path.join(OUT, "files")
U64 = np.uint64


# ---- generated arrays -------------------------------------------------------------------------

def splitmix(seed, n):
    i = np.arange(1, n + 1, dtype=np.uint64)
    with np.errstate(over="ignore"):
        z = U64(seed) + i * U64(0x9E3779B97F4A7C15)
        z = (z ^ (z >> U64(30))) * U64(0xBF58476D1CE4E5B9)
        z = (z ^ (z >> U64(27))) * U64(0x94D049BB133111EB)
        z = z ^ (z >> U64(31))
    return z


def generate(spec):
    n = spec["n"]
    kind = spec["kind"]
    if kind == "uniform":
        u = (splitmix(spec["seed"], n) >> U64(40)).astype(np.float64) / 2.0 ** 24
        return ((u - 0.5) * spec["scale"] + spec["offset"]).astype(np.float32)
    if kind == "bits":
        z = splitmix(spec["seed"], n)
        span = spec["emax"] - spec["emin"] + 1
        e = U64(spec["emin"]) + (z >> U64(32)) % U64(span)
        bits = (z & U64(0x807FFFFF)) | (e << U64(23))
        return bits.astype(np.uint32).view(np.float32)
    if kind == "constant":
        return np.full(n, spec["value"], dtype=np.float32)
    if kind == "ramp":
        return (spec["start"] + np.arange(n, dtype=np.float64) * spec["step"]).astype(np.float32)
    raise ValueError(kind)


def f32bits(x):
    return int(np.float32(x).view(np.uint32))


def sha(data):
    return hashlib.sha256(np.ascontiguousarray(data, dtype="<f4").tobytes()).hexdigest()


def floats(values):
    return [float(v) for v in np.asarray(values, dtype=np.float64).ravel()]


def stats_bits(dmap):
    d = dmap.data
    return [f32bits(d.min()), f32bits(d.max()), f32bits(d.mean()), f32bits(d.std())]


# ---- map files --------------------------------------------------------------------------------

CODES = {0: "i1", 1: "i2", 2: "f4", 6: "u2", 12: "f2"}


def header(nc, nr, ns, mode=2, endian="<", nstart=(0, 0, 0), mxyz=None, cell=None, voxel=(1.0, 1.0, 1.0),
           angles=(90.0, 90.0, 90.0), mapcrs=(1, 2, 3), nsymbt=0, origin=(0.0, 0.0, 0.0), stamp=None,
           tag=b"MAP ", imod=None):
    """An MRC header in the manner of MDATK's tests/test_density.py::_raw_map, every field settable."""
    ints = np.zeros(256, dtype=f"{endian}i4")
    fl = ints.view(f"{endian}f4")
    ints[0:4] = (nc, nr, ns, mode)
    ints[4:7] = nstart
    n_cart = [0, 0, 0]
    for k, m in enumerate(mapcrs):
        if 1 <= m <= 3:
            n_cart[m - 1] = (nc, nr, ns)[k]
    ints[7:10] = n_cart if mxyz is None else mxyz
    fl[10:13] = np.array(n_cart) * np.array(voxel) if cell is None else cell
    fl[13:16] = angles
    ints[16:19] = mapcrs
    ints[23] = nsymbt
    ints[27] = 20140
    fl[49:52] = origin
    if imod is not None:
        ints[38], ints[39] = imod
    raw = bytearray(ints.tobytes())
    raw[208:212] = tag
    raw[212:214] = stamp if stamp is not None else (b"\x44\x44" if endian == "<" else b"\x11\x11")
    return bytes(raw)


def body(values, mode=2, endian="<", unsigned=False):
    code = "u1" if (mode == 0 and unsigned) else CODES[mode]
    return np.ascontiguousarray(values).astype(f"{endian}{code}").tobytes()


def mrc(shape_crs=(5, 4, 3), seed=1, mode=2, endian="<", unsigned=False, values=None, extra=b"", **kw):
    """Header + extended header + data; shape_crs = (nc, nr, ns)."""
    nc, nr, ns = shape_crs
    n = nc * nr * ns
    if values is None:
        if mode == 2:
            values = generate({"kind": "uniform", "seed": seed, "n": n, "scale": 4.0, "offset": 0.5})
        elif mode == 12:
            values = generate({"kind": "uniform", "seed": seed, "n": n, "scale": 100.0, "offset": 0.0})
        else:
            lo, hi = {0: (0, 255) if unsigned else (-128, 127), 1: (-32768, 32767), 6: (0, 65535)}[mode]
            values = lo + (splitmix(seed, n) % U64(hi - lo + 1)).astype(np.int64)
    head = header(nc, nr, ns, mode=mode, endian=endian, **kw)
    nsymbt = kw.get("nsymbt", 0)
    return head + bytes(range(256)) * (max(nsymbt, 0) // 256) + b"\x5a" * (max(nsymbt, 0) % 256) + \
        body(values, mode, endian, unsigned) + extra


def gz(data):
    return gzip.compress(data, mtime=0)


def raw_deflate(data):
    c = zlib.compressobj(9, zlib.DEFLATED, -15)
    return c.compress(data) + c.flush()


def member(data, flags=0, extra=b"", name=b"", comment=b"", crc=None, isize=None, method=8):
    head = bytes([0x1F, 0x8B, method, flags]) + b"\0\0\0\0" + b"\x02\x03"
    if flags & 4:
        head += struct.pack("<H", len(extra)) + extra
    if flags & 8:
        head += name + b"\0"
    if flags & 16:
        head += comment + b"\0"
    if flags & 2:
        head += b"\xab\xcd"  # header CRC: not checked by Python
    crc = zlib.crc32(data) if crc is None else crc
    isize = len(data) & 0xFFFFFFFF if isize is None else isize
    return head + raw_deflate(data) + struct.pack("<II", crc, isize)


def stored_member(data):
    """gzip member of stored (uncompressed) deflate blocks."""
    out = bytes([0x1F, 0x8B, 8, 0]) + b"\0" * 4 + b"\0\x03"
    pos = 0
    while True:
        chunk = data[pos:pos + 700]
        pos += len(chunk)
        last = pos >= len(data)
        out += bytes([1 if last else 0]) + struct.pack("<HH", len(chunk), len(chunk) ^ 0xFFFF) + chunk
        if last:
            break
    return out + struct.pack("<II", zlib.crc32(data), len(data))


def build_files():
    files = []  # (name, bytes, deviation)

    def add(name, data, deviation=None):
        files.append((name, data, deviation))

    # every axis order in both byte orders, with an extended header and an ORIGIN
    for mapcrs in ((1, 2, 3), (1, 3, 2), (2, 1, 3), (2, 3, 1), (3, 1, 2), (3, 2, 1)):
        for endian in "<>":
            tag = "".join(map(str, mapcrs)) + ("le" if endian == "<" else "be")
            add(f"axes_{tag}.mrc", mrc(seed=sum(mapcrs) * 7 + len(files), endian=endian, mapcrs=mapcrs,
                                       nsymbt=64, voxel=(0.8, 1.1, 1.3), origin=(10.0, -5.0, 2.0)))
            add(f"axes_nstart_{tag}.mrc", mrc(seed=len(files), endian=endian, mapcrs=mapcrs,
                                              nstart=(-4, 2, 10), voxel=(2.0, 1.5, 0.5)))
    # modes
    for endian in "<>":
        e = "le" if endian == "<" else "be"
        add(f"mode0_signed_{e}.mrc", mrc(mode=0, endian=endian, seed=3))
        add(f"mode0_notag_{e}.mrc", mrc(mode=0, endian=endian, seed=3, unsigned=True, tag=b"CCP4"))
        add(f"mode1_{e}.mrc", mrc(mode=1, endian=endian, seed=4,
                                  values=np.array([-32768, 32767, 0, -1, 1] + list(range(-100, 1100, 22)))))
        add(f"mode6_{e}.mrc", mrc(mode=6, endian=endian, seed=5,
                                  values=np.array([0, 65535, 1, 32768] + list(range(0, 56000, 1000)))))
        halves = np.array([0.0, -0.0, 65504.0, -65504.0, 6.1e-5, 5.96e-8, -5.96e-8, 1e-6, 0.1, -2.5, 1000.0,
                           np.inf, -np.inf] + list(np.linspace(-3, 3, 47)), dtype=np.float16)
        add(f"mode12_{e}.mrc", mrc(mode=12, endian=endian, values=halves))
        add(f"mode2_specials_{e}.mrc", mrc(mode=2, endian=endian, values=np.array(
            [0.0, -0.0, 1e-45, -1e-45, 1.17549435e-38, 3.4028235e38, -3.4028235e38, 1.5, -2.25]
            + [0.125 * k for k in range(51)], dtype=np.float32)))
    add("mode0_imod_unsigned.mrc", mrc(mode=0, seed=6, unsigned=True, imod=(1146047817, 0)))
    add("mode0_imod_flag1.mrc", mrc(mode=0, seed=6, imod=(1146047817, 1)))
    add("mode0_imod_flag3.mrc", mrc(mode=0, seed=6, imod=(1146047817, 3)))
    add("mode0_imod_flag_neg2.mrc", mrc(mode=0, seed=6, unsigned=True, imod=(1146047817, -2)))
    add("mode0_imod_other_stamp.mrc", mrc(mode=0, seed=6, imod=(1146047818, 0)))
    add("mode0_imod_notag.mrc", mrc(mode=0, seed=6, unsigned=True, tag=b"\0\0\0\0", imod=(1146047817, 1)))
    add("mode2_nan.mrc", mrc(values=np.array([1.0, np.nan] + [2.0] * 58, dtype=np.float32)))
    add("mode2_inf.mrc", mrc(values=np.array([1.0, np.inf] + [2.0] * 58, dtype=np.float32)))
    add("mode2_constant.mrc", mrc(values=np.full(60, 0.1, dtype=np.float32)))
    for mode in (3, 4, 101, -1):
        add(f"mode_unsupported_{mode}.mrc", header(5, 4, 3, mode=mode) + b"\0" * 240)
    # machine stamps and byte-order guessing
    add("stamp_4441.mrc", mrc(seed=8, stamp=b"\x44\x41"))
    add("stamp_1111_be.mrc", mrc(seed=8, endian=">", stamp=b"\x11\x11"))
    add("stamp_none_le.mrc", mrc(seed=8, stamp=b"\0\0"))
    add("stamp_none_be.mrc", mrc(seed=8, endian=">", stamp=b"\0\0"))
    add("stamp_other_be_mode0.mrc", mrc(mode=0, seed=8, endian=">", stamp=b"\x41\x41"))
    add("stamp_none_be_mode1_crs0.mrc", mrc(mode=1, seed=8, endian=">", stamp=b"\0\0", mapcrs=(0, 0, 0)))
    add("stamp_none_mode0_crs0.mrc", mrc(mode=0, seed=8, endian=">", stamp=b"\0\0", mapcrs=(0, 0, 0)))
    add("stamp_wrong_le_for_be.mrc", mrc(seed=8, endian=">", stamp=b"\x44\x44"))
    add("stamp_none_invalid.mrc", header(5, 4, 3, mode=3, stamp=b"\0\0") + b"\0" * 240)
    add("stamp_none_bad_crs.mrc", header(5, 4, 3, stamp=b"\x20\0", mapcrs=(1, 2, 4)) + b"\0" * 240)
    # NSTART / ORIGIN / MXYZ / cell / angles
    add("origin_zero_nstart.mrc", mrc(nstart=(-4, 2, 10), voxel=(2.0, 2.0, 2.0)))
    add("origin_nan.mrc", mrc(nstart=(1, 2, 3), origin=(np.nan, 1.0, 2.0)))
    add("origin_inf.mrc", mrc(nstart=(1, 2, 3), origin=(np.inf, 0.0, 0.0)))
    add("origin_negzero.mrc", mrc(nstart=(1, 2, 3), origin=(-0.0, 0.0, 0.0)))
    add("origin_tiny.mrc", mrc(nstart=(1, 2, 3), origin=(0.0, 0.0, 1e-30)))
    add("origin_large.mrc", mrc(origin=(1e7, -3.3e5, 12.345678)))
    add("cell_zero.mrc", mrc(nstart=(3, -3, 1), cell=(0.0, 10.0, 10.0)))
    add("cell_negative.mrc", mrc(nstart=(3, -3, 1), cell=(5.0, -10.0, 10.0)))
    add("cell_nan.mrc", mrc(cell=(np.nan, 10.0, 10.0)))
    add("mxyz_zero.mrc", mrc(mxyz=(0, 0, 0), cell=(10.0, 8.0, 6.0)))
    add("mxyz_negative.mrc", mrc(mxyz=(-1, 4, -7), cell=(10.0, 8.0, 6.0)))
    add("mxyz_other.mrc", mrc(mxyz=(10, 20, 30), cell=(10.0, 8.0, 6.0), nstart=(1, 1, 1)))
    add("mxyz_permuted.mrc", mrc(mapcrs=(3, 1, 2), mxyz=(0, 9, 0), cell=(10.0, 8.0, 6.0), nstart=(5, 6, 7)))
    add("crs_000.mrc", mrc(mapcrs=(0, 0, 0), nstart=(1, 2, 3)))
    add("angles_120.mrc", mrc(angles=(90.0, 90.0, 120.0)))
    add("angles_slight.mrc", mrc(angles=(90.02, 90.0, 90.0)))
    add("angles_within.mrc", mrc(angles=(90.005, 90.0, 89.995)))
    add("angles_zero.mrc", mrc(angles=(0.0, 0.0, 0.0)))
    add("angles_negative.mrc", mrc(angles=(-90.0, 120.0, 90.0)))
    add("angles_nan.mrc", mrc(angles=(np.nan, 90.0, 90.0)))
    add("angles_60.mrc", mrc(angles=(60.0, 60.0, 60.0)))
    # dimensions and axis orders that are rejected
    add("dims_zero.mrc", header(0, 4, 3) + b"\0" * 240)
    add("dims_negative.mrc", header(5, -1, 3, mapcrs=(1, 2, 3)) + b"\0" * 240)
    for crs in ((1, 1, 3), (1, 2, 4), (0, 0, 1), (3, 2, 0), (-1, 2, 3)):
        add("crs_bad_" + "_".join(map(str, crs)).replace("-", "m") + ".mrc",
            header(5, 4, 3, mapcrs=crs) + b"\0" * 240)
    # extended headers, lengths
    add("nsymbt_1000.mrc", mrc(nsymbt=1000, seed=9))
    add("nsymbt_negative.mrc", mrc(nsymbt=-5, seed=9))
    add("nsymbt_beyond.mrc", header(5, 4, 3, nsymbt=2 ** 31 - 1) + b"\0" * 300)
    add("trailing_bytes.mrc", mrc(seed=9, extra=b"extra bytes after the data"))
    full = mrc(seed=10)
    add("short_by_one.mrc", full[:-1])
    add("header_only.mrc", full[:1024])
    add("short_header.mrc", full[:1023])
    add("empty.mrc", b"")
    add("one_byte.mrc", b"\x1f")
    add("zeros_100.mrc", b"\0" * 100)
    add("zeros_2000.mrc", b"\0" * 2000)
    add("huge_dims.mrc", header(2 ** 31 - 1, 2 ** 31 - 1, 2 ** 31 - 1) + b"\0" * 64, "shorter")
    add("huge_dims_mode0.mrc", header(100000, 100000, 1000, mode=0) + b"\0" * 64, "shorter")
    # gzip
    m = mrc(seed=11, nstart=(-2, 1, 4), voxel=(1.5, 1.5, 1.5))
    add("gzip.mrc.gz", gz(m))
    add("gzip_members.map.gz", gz(m[:500]) + gz(m[500:1100]) + gz(m[1100:]))
    add("gzip_empty_member.map.gz", gz(b"") + gz(m))
    add("gzip_flags.map.gz", member(m[:900], flags=4 | 8 | 16 | 2, extra=b"xyz\0", name=b"map.mrc",
                                    comment=b"a comment") + member(m[900:], flags=8, name=b""))
    add("gzip_padding.map.gz", gz(m[:1500]) + b"\0" * 17 + gz(m[1500:]))
    add("gzip_stored.map.gz", stored_member(m))
    add("gzip_trailing_garbage.map.gz", gz(m) + b"garbage")
    add("gzip_trailer_cut.map.gz", gz(m)[:-4])
    add("gzip_magic_only.map.gz", b"\x1f\x8b", "eof")
    add("gzip_cut_header.map.gz", gz(m)[:200], "eof")
    add("gzip_cut_data.map.gz", gz(m)[:-40], "eof")
    add("gzip_cut_member.map.gz", gz(m[:1100]) + gz(m[1100:])[:12], "eof")
    bad = bytearray(gz(m))
    bad[10] = 0xFF  # first deflate block: type 3
    add("gzip_bad_block.map.gz", bytes(bad), "zlib")
    add("gzip_bad_stored.map.gz", stored_member(m)[:11] + b"\x10\x00\x10\x00" + stored_member(m)[15:], "zlib")
    add("gzip_bad_crc.map.gz", member(m[:1200], crc=0x12345678) + gz(m[1200:]), "gzip")
    add("gzip_bad_isize.map.gz", member(m[:1200], isize=7) + gz(m[1200:]), "gzip")
    add("gzip_bad_magic.map.gz", gz(m[:1200]) + b"XY" + gz(m[1200:]), "gzip")
    add("gzip_bad_magic_quote.map.gz", gz(m[:1200]) + b"'x" + gz(m[1200:]), "gzip")
    add("gzip_one_byte_tail.map.gz", gz(m[:1200]) + b"\x8b", "gzip")
    add("gzip_bad_method.map.gz", gz(m[:1200]) + member(m[1200:], method=7), "gzip")
    add("gzip_huge_dims.map.gz", gz(header(2 ** 31 - 1, 2 ** 31 - 1, 2 ** 31 - 1)), "shorter")
    add("gzip_short.mrc.gz", gz(m[:1000]))
    # a larger map (123 kB, a pattern of period 101 that compresses well): chunk boundaries of the
    # decompression
    q = (np.arange(40 * 32 * 24) * 7919 % 101) * 0.25 - 3.0
    big = mrc((40, 32, 24), values=q.astype(np.float32), mapcrs=(2, 3, 1), voxel=(1.2, 1.2, 1.2))
    add("large.map.gz", gz(big))
    add("large_members.map.gz", b"".join(gz(big[k:k + 9000]) for k in range(0, len(big), 9000)))
    return files


@contextlib.contextmanager
def memory_limit(nbytes):
    soft, hard = resource.getrlimit(resource.RLIMIT_AS)
    resource.setrlimit(resource.RLIMIT_AS, (nbytes, hard))
    try:
        yield
    finally:
        resource.setrlimit(resource.RLIMIT_AS, (soft, hard))


def available(data):
    """Bytes Python's gzip reader returns from truncated data before its EOFError."""
    reader = gzip.GzipFile(fileobj=io.BytesIO(data))
    count = 0
    try:
        while reader.read(1):
            count += 1
    except EOFError:
        pass
    return count


def read_cases(files):
    if os.path.isdir(FILES):
        shutil.rmtree(FILES)
    os.makedirs(FILES)
    cases = []
    for name, data, deviation in files:
        path = os.path.join(FILES, name)
        with open(path, "wb") as fh:
            fh.write(data)
        case = {"name": name}
        try:
            with memory_limit(4 << 30):
                dmap = read_map(path)
            if deviation is not None:
                raise SystemExit(f"{name}: expected an error ({deviation})")
            case.update({"shape_xyz": list(dmap.shape_xyz), "origin": floats(dmap.origin),
                         "voxel": floats(dmap.voxel), "sha256": sha(dmap.data),
                         "stats": stats_bits(dmap), "first": floats(dmap.data.ravel()[:4])})
        except (MapError, EOFError, OSError, zlib.error, MemoryError, OverflowError) as exc:
            message = str(exc).replace(path, "{path}")
            if isinstance(exc, MapError):
                if deviation not in (None, "shorter"):
                    raise SystemExit(f"{name}: MapError {message!r}, expected {deviation}")
                case["error"] = message
            else:
                case["python"] = {"type": type(exc).__name__, "message": message}
                if deviation == "eof":  # truncated gzip: where the data ends decides the message
                    deviation = "header" if available(data) < 1024 else "shorter"
                expected = {"header": "{path}: file too short for an MRC header",
                            "shorter": "{path}: file is shorter than its header says",
                            "zlib": "{path}: " + message, "gzip": "{path}: " + message}.get(deviation)
                if expected is None:
                    raise SystemExit(f"{name}: unexpected {type(exc).__name__}: {message}")
                kinds = {"header": (EOFError,), "shorter": (EOFError, MemoryError, OverflowError),
                         "zlib": (zlib.error,), "gzip": (gzip.BadGzipFile,)}[deviation]
                if not isinstance(exc, kinds):
                    raise SystemExit(f"{name}: {type(exc).__name__} is not a {deviation} case")
                case["error"] = expected
        cases.append(case)
    return cases


def info_cases(names):
    cases = []
    for name in names:
        path = os.path.join(FILES, name)
        out = io.StringIO()
        with contextlib.redirect_stdout(out):
            assert cli.main(["density", "info", path]) == 0
        cases.append({"name": name, "output": out.getvalue().replace(os.path.abspath(path), "{path}")})
    return cases


# ---- statistics -------------------------------------------------------------------------------

def stats_cases():
    specs = []
    sizes = list(range(1, 21)) + [63, 64, 65, 127, 128, 129, 135, 136, 137, 255, 256, 257, 1000, 1023,
                                  1024, 1025, 4096, 8191, 8192, 8193, 65535, 65536, 65537, 999_983]
    for k, n in enumerate(sizes):
        specs.append({"kind": "uniform", "seed": 100 + k, "n": n, "scale": 10.0 ** (k % 7 - 2),
                      "offset": [0.0, 1.0, -3.5, 1000.0][k % 4]})
    for k, n in enumerate([1, 7, 8, 9, 100, 129, 1000, 4097, 65537, 300_001]):
        specs.append({"kind": "bits", "seed": 200 + k, "n": n, "emin": 100, "emax": 154})
        specs.append({"kind": "bits", "seed": 300 + k, "n": n, "emin": 1, "emax": 254})
    for n, v in ((1, 0.1), (3, -0.0), (17, 0.1), (1000, 0.1), (100_000, 0.1), (16_777_217, 0.1),
                 (5, 3.4e38), (300, 1e30)):
        specs.append({"kind": "constant", "n": n, "value": v})
    for n, start, step in ((10, 0.0, 1.0), (1000, -500.0, 1.0), (100_000, 0.0, 0.001),
                           (3 * 2 ** 22 + 5, -1.0, 1e-7), (16_777_217, 0.0, 1.0)):
        specs.append({"kind": "ramp", "n": n, "start": start, "step": step})
    specs.append({"kind": "uniform", "seed": 999, "n": 2 ** 24 + 1, "scale": 3.0, "offset": 0.25})
    cases = []
    for spec in specs:
        dmap = DensityMap(generate(spec).reshape(1, 1, -1), np.zeros(3), np.ones(3))
        levels = [(s, D.map_level(dmap, s)) for s in (0.0, 1.0, 2.0, -1.5, 3.7)]
        cases.append({"data": spec, "stats": stats_bits(dmap),
                      "levels": [[s, float(v)] for s, v in levels]})
    return cases


# ---- writing ----------------------------------------------------------------------------------

def write_cases():
    cases = []
    tmp = os.path.join(OUT, "_write.mrc")
    shapes = [(1, 1, 1), (7, 6, 5), (16, 12, 9), (64, 48, 40), (33, 1, 2)]
    origins = [(1.5, -2.0, 3.25), (-3.0, 1.5, 6.0), (0.0, 0.0, 0.0), (1.10003, 2.2, -3.3), (1e12, 5.0, 0.0),
               (12.000049, -7.99995, 0.0), (np.inf, 0.0, 0.0), (-0.0, 0.5, -0.0)]
    voxels = [(1.1, 1.1, 1.1), (1.5, 1.5, 1.5), (0.8, 1.2, 3.0), (1.1, 1.1, 1.1), (1.0, 1.0, 1.0),
              (2.0, 2.0, 2.0), (1.0, 1.0, 1.0), (0.25, 0.5, 1.0)]
    for k, (origin, voxel) in enumerate(zip(origins, voxels)):
        for j, shape in enumerate(shapes):
            if (k + j) % 2 and shape != (64, 48, 40):
                continue
            nx, ny, nz = shape
            spec = {"kind": "uniform", "seed": 500 + 10 * k + j, "n": nx * ny * nz,
                    "scale": 2.0 + k, "offset": 0.5 * j}
            if k == 3 and j == 2:
                spec = {"kind": "bits", "seed": 77, "n": nx * ny * nz, "emin": 90, "emax": 160}
            data = generate(spec).reshape(nz, ny, nx)
            write_map(tmp, DensityMap(data, np.array(origin, dtype=float), np.array(voxel, dtype=float)))
            with open(tmp, "rb") as fh:
                raw = fh.read()
            assert raw[224:1024] == b"\0" * 800
            back = read_map(tmp)
            cases.append({"data": spec, "shape_xyz": list(shape), "origin": floats(origin),
                          "voxel": floats(voxel), "sha256": hashlib.sha256(raw).hexdigest(),
                          "size": len(raw), "header": raw[:224].hex(),
                          "read": {"origin": floats(back.origin), "voxel": floats(back.voxel),
                                   "sha256": sha(back.data)}})
    os.remove(tmp)
    try:
        write_map(tmp, DensityMap(np.zeros((0, 2, 2), np.float32), np.zeros(3), np.ones(3)))
        raise SystemExit("empty map written")
    except ValueError as exc:
        empty_error = str(exc)
    return {"cases": cases, "empty_error": empty_error, "empty_written": os.path.exists(tmp)}


# ---- simulated maps ---------------------------------------------------------------------------

class LibmExp:
    """numpy with exp computed by the C library (math.exp), element by element."""

    def __getattr__(self, name):
        return getattr(np, name)

    @staticmethod
    def exp(x):
        x = np.asarray(x, dtype=np.float64)
        return np.array([math.exp(v) for v in x.ravel()], dtype=np.float64).reshape(x.shape)


def simulate_cases():
    rng = np.random.default_rng(7)
    elements = ["C", " n ", "o", "Mg", "NA", "Zn", "X", "", "cl", "ǆ", "ß", "P\t", "fe", "se", "BA",
                "H", "s", "Cs"]
    weights_case = {"elements": elements, "weights": floats(D.atom_weights(elements))}
    specs = [
        {"xyz": [[0.0, 0.0, 0.0], [12.0, 0.0, 0.0]], "weights": [6.0, 6.0], "resolution": 4.0, "spacing": 0.5},
        {"xyz": [[0.3, 0.7, -1.1]], "weights": [6.0], "resolution": 3.0},
        {"xyz": [[0.3, 0.7, -1.1]], "weights": [6.0], "resolution": 3.0, "spacing": 0.0},
        {"xyz": [[1.0, 2.0, 3.0], [1.5, 2.5, 3.5]], "weights": [8.0, 1.0], "resolution": 2.0,
         "sigma_factor": 0.3, "pad_sigmas": 3.0},
        {"xyz": [[-1000.25, 2000.5, -3.75]], "weights": [15.0], "resolution": 5.0, "spacing": 1.3},
        {"xyz": [[0.0, 0.0, 0.0]], "weights": [1.0], "resolution": 0.3, "spacing": 1.0},
    ]
    for k, (n, res, spacing) in enumerate(((40, 6.0, 1.5), (40, 8.0, 1.6), (25, 1.0, None), (300, 4.0, None))):
        xyz = rng.normal(size=(n, 3)) * 6 + rng.uniform(-50, 50, size=3)
        els = [elements[i] for i in rng.integers(0, len(elements), size=n)]
        spec = {"xyz": [floats(p) for p in xyz], "elements": els, "resolution": res}
        if spacing is not None:
            spec["spacing"] = spacing
        specs.append(spec)
    cases = []
    for spec in specs:
        xyz = np.array(spec["xyz"], dtype=np.float64)
        w = D.atom_weights(spec["elements"]) if "elements" in spec else np.array(spec["weights"])
        kwargs = {k: spec[k] for k in ("spacing", "sigma_factor", "pad_sigmas") if k in spec}
        real = D.simulate_map(xyz, w, spec["resolution"], **kwargs)
        D.np = LibmExp()
        try:
            sim = D.simulate_map(xyz, w, spec["resolution"], **kwargs)
        finally:
            D.np = np
        case = dict(spec)
        case.update({"shape_xyz": list(sim.shape_xyz), "origin": floats(sim.origin), "voxel": floats(sim.voxel),
                     "sha256": sha(sim.data), "stats": stats_bits(sim),
                     "numpy_sha256": sha(real.data),
                     "numpy_differing_voxels": int(np.count_nonzero(real.data != sim.data))})
        cases.append(case)
    errors = []
    for xyz, res in ((np.zeros((0, 3)), 3.0), (np.zeros((1, 3)), 0.0), (np.zeros((1, 3)), -1.0)):
        try:
            D.simulate_map(xyz, np.ones(len(xyz)), res)
        except D.FitError as exc:
            errors.append({"atoms": len(xyz), "resolution": res, "error": str(exc)})
    return {"sigma_factor": D.SIGMA_FACTOR, "atom_weights": weights_case, "cases": cases, "errors": errors}


def dump(path, data):
    """Lists of cases one per line."""
    def render(value, indent):
        if isinstance(value, list) and value and isinstance(value[0], dict):
            rows = ",\n".join(" " * (indent + 1) + json.dumps(v, separators=(",", ":")) for v in value)
            return "[\n" + rows + "\n" + " " * indent + "]"
        if isinstance(value, dict):
            items = ",\n".join(" " * (indent + 1) + json.dumps(k) + ": " + render(v, indent + 1)
                               for k, v in value.items())
            return "{\n" + items + "\n" + " " * indent + "}"
        return json.dumps(value, separators=(",", ":"))
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(render(data, 0) + "\n")


def main():
    os.makedirs(OUT, exist_ok=True)
    files = build_files()
    reads = read_cases(files)
    dump(os.path.join(OUT, "read_cases.json"), {"cases": reads})
    ok = [c["name"] for c in reads if "error" not in c]
    dump(os.path.join(OUT, "info_cases.json"),
         {"cases": info_cases([n for n in ok if n.startswith(("axes_123le", "origin_", "mode12", "large",
                                                              "cell_zero", "gzip.mrc"))])})
    dump(os.path.join(OUT, "stats_cases.json"), {"cases": stats_cases()})
    dump(os.path.join(OUT, "write_cases.json"), write_cases())
    dump(os.path.join(OUT, "simulate_cases.json"), simulate_cases())
    print(len(reads), "read cases,", sum("error" in c for c in reads), "errors;",
          sum(len(d) for _n, d, _e in files), "bytes of files")


main()
