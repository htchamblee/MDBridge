"""What MDAnalysis 2.11.0 reads from coordinate files for MDATK (tests/data/io/trajectories).

    /opt/mdbridge-venv/bin/python -I -B make_trajectory_cases.py MDATK_CHECKOUT OUT_DIR

Needs MDAnalysis 2.11.0, numpy and scipy (MDAnalysis' NetCDF backend without netCDF4); MDATK's
checkout provides mdatk.demo and the extension table of mdatk.io.loader (nothing is written into
it). Writes

  files/*               small trajectories: NetCDF written by MDAnalysis' NCDFWriter and crafted
                        with scipy (CDF-1 / CDF-2, every variable type, scale factors, no time /
                        no box, malformed Amber conventions), Amber mdcrd and rst7 texts, DCD files
                        (MDAnalysis' writer and crafted: endianness, X-PLOR, unit cells as cosines /
                        degrees / box vectors, 4D, fixed atoms, truncation), XTC / TRR (MDAnalysis'
                        writers: compressed and uncompressed, precisions, wide coordinate ranges,
                        triclinic boxes; crafted double-precision TRR) and broken files
  demo_dry.nc           mdatk.demo.build_demo(n_frames=24, solvent=False)'s trajectory (the
                        tests' demo_dry fixture; tests/data/topology/demo_dry.parm7 is its topology)
  writer/*.nc           NCDFWriter's output for frames of the cases (with and without boxes; the
                        demo trajectory rewritten is demo_dry.nc itself)
  cases.json            per case: the files (chained when several), format, dt and atom count given
                        to the Universe, and MDAnalysis' n_frames, dt, every frame's time, box and
                        positions (float32 bits; SHA-256 of the little-endian float32 bytes for the
                        demo), the TRJ header; or the error {"type", "message"} with each file's
                        path written {path0}, {path1}, ...
  formats.json          MDATK's extension table, MDAnalysis' per-file format guess and its errors

Times and dt are float.hex() strings, boxes and positions big-endian float32 bits as hex.
"""

import hashlib
import json
import math
import os
import shutil
import struct
import sys
import tempfile
import warnings

warnings.filterwarnings("ignore")

if len(sys.argv) != 3:
    sys.exit("usage: make_trajectory_cases.py MDATK_CHECKOUT OUT_DIR")
sys.path.insert(0, os.path.abspath(sys.argv[1]))

import MDAnalysis as mda  # noqa: E402
import numpy as np  # noqa: E402
import scipy.io  # noqa: E402
from MDAnalysis.coordinates.memory import MemoryReader  # noqa: E402

from mdatk.demo import build_demo  # noqa: E402
from mdatk.io import loader  # noqa: E402

OUT = os.path.abspath(sys.argv[2])
FILES = os.path.join(OUT, "files")
WRITER = os.path.join(OUT, "writer")
TMP = tempfile.mkdtemp(prefix="trajectory-cases-")
for d in (FILES, WRITER):
    shutil.rmtree(d, ignore_errors=True)
    os.makedirs(d)

N = 7
rng = np.random.default_rng(2026)
COORDS = (rng.random((4, N, 3)) * 40 - 20).astype(np.float32)
BOX = np.array([30.5, 31.25, 32.125, 90, 90, 90], dtype=np.float32)
TRICLINIC = np.array([30.5, 31.25, 32.125, 80, 95.5, 110.25], dtype=np.float32)


def f32hex(x):
    return struct.pack(">f", float(x)).hex()


def tmp(name):
    return os.path.join(TMP, name)


def universe(coords, dt=1.0, box=None, times=None):
    u = mda.Universe.empty(coords.shape[1], trajectory=True)
    u.load_new(coords.copy(), format=MemoryReader, dt=dt, dimensions=box)
    return u


def mda_write(name, coords, dt=1.0, box=None, fmt=None, **kw):
    u = universe(coords, dt, box)
    with mda.Writer(tmp(name), n_atoms=coords.shape[1], format=fmt, **kw) as w:
        for _ in u.trajectory:
            w.write(u.atoms)
    return name


def raw(name, data):
    with open(tmp(name), "wb") as f:
        f.write(data)
    return name


# ---------------------------------------------------------------- NetCDF

def nc_crafted(name, *, version=2, coords=COORDS, conventions="AMBER", conv_version="1.0", time="f",
               times=None, box=True, coord_type="f", coord_scale=None, time_scale=None,
               lengths_scale=None, angles_scale=None, units=None, extras=False, drop=(),
               spatial=3, natom_dim=None, gattrs=None, attrs=None, frame_fixed=False, angles=True):
    f = scipy.io.netcdf_file(tmp(name), "w", version=version)
    if conventions is not None:
        f.Conventions = conventions
    if conv_version is not None:
        f.ConventionVersion = conv_version
    f.program = "make_trajectory_cases"
    f.programVersion = "1"
    for k, v in (gattrs or {}).items():
        setattr(f, k, v)
    nf, na = coords.shape[:2]
    if "frame" not in drop:
        f.createDimension("frame", nf if frame_fixed else None)
    if "atom" not in drop:
        f.createDimension("atom", natom_dim if natom_dim is not None else na)
    if "spatial" not in drop:
        f.createDimension("spatial", spatial)
    f.createDimension("cell_spatial", 3)
    f.createDimension("cell_angular", 3)
    f.createDimension("label", 5)
    units = dict(units or {})
    if extras:
        f.createDimension("two", 2)
        for t in "bhifd":
            v = f.createVariable("fixed_" + t, t, ("atom", "two"))
            v[:] = (np.arange(na * 2).reshape(na, 2) % 100).astype(t)
            v.note = "x" * (1 + "bhifd".index(t))
            r = f.createVariable("rec_" + t, t, ("frame",))
            r[:] = np.arange(nf).astype(t)
        c = f.createVariable("fixed_c", "c", ("two",))
        c[:] = np.array([b"a", b"b"])
        rc = f.createVariable("rec_c", "c", ("frame", "label"))
        rc[:] = np.array([list(b"hello")] * nf, dtype="S1")
    if "spatial" not in drop:
        sp = f.createVariable("spatial", "c", ("spatial",))
        if spatial == 3:
            sp[:] = np.array(list("xyz"), dtype="S1")
    if "coordinates" not in drop:
        dims = tuple(d for d in ("frame", "atom", "spatial") if d not in drop)
        v = f.createVariable("coordinates", coord_type, dims)
        if "coordinates.units" not in drop:
            v.units = units.get("coordinates", "angstrom")
        data = coords.astype(np.float64)
        if coord_scale is not None:
            v.scale_factor = coord_scale
            data = data / float(coord_scale)
        if coord_type in "bhi":
            data = np.round(data)
        if "atom" in drop:
            data = data[:, 0, :]
        if "frame" in drop:
            data = coords[0].astype(np.float64)
        v[:] = data.astype(coord_type)
        for k, val in (attrs or {}).get("coordinates", {}).items():
            setattr(v, k, val)
    if time:
        t = f.createVariable("time", time, ("frame",))
        t.units = units.get("time", "picosecond")
        values = np.asarray(times if times is not None else np.arange(nf) * 0.75 + 2.0, dtype=np.float64)
        if time_scale is not None:
            t.scale_factor = time_scale
        t[:] = values.astype(time)
    if box:
        lv = f.createVariable("cell_lengths", "d", ("frame", "cell_spatial"))
        lv.units = units.get("cell_lengths", "angstrom")
        lv[:] = np.array([[30.1 + k, 31.7, 32.3] for k in range(nf)]).reshape(nf, 3)
        if lengths_scale is not None:
            lv.scale_factor = lengths_scale
        if angles:
            av = f.createVariable("cell_angles", "d", ("frame", "cell_angular"))
            av.units = units.get("cell_angles", "degree")
            av[:] = np.array([[90.0, 91.3, 119.9]] * nf).reshape(nf, 3)
            if angles_scale is not None:
                av.scale_factor = angles_scale
    for var, values in (attrs or {}).items():
        if var != "coordinates" and var in f.variables:
            for k, val in values.items():
                setattr(f.variables[var], k, val)
    f.close()
    return name


def netcdf_files():
    out = {}
    out["mda_box"] = mda_write("mda_box.nc", COORDS, dt=0.5, box=BOX)
    out["mda_nobox"] = mda_write("mda_nobox.nc", COORDS, dt=0.1)
    out["mda_triclinic"] = mda_write("mda_triclinic.nc", COORDS, dt=2.0, box=TRICLINIC)
    out["mda_one_frame"] = mda_write("mda_one_frame.nc", COORDS[:1], dt=0.5, box=BOX)
    out["mda_scaled"] = mda_write("mda_scaled.nc", COORDS, dt=0.5, box=TRICLINIC, scale_time=2.0,
                                  scale_coordinates=0.1, scale_cell_lengths=0.5, scale_cell_angles=1.5)
    out["cdf1"] = nc_crafted("cdf1.nc", version=1)
    out["cdf2"] = nc_crafted("cdf2.nc")
    out["no_time"] = nc_crafted("no_time.nc", time=None)
    out["no_time_no_box"] = nc_crafted("no_time_no_box.nc", time=None, box=False)
    out["only_coordinates_byte"] = nc_crafted("only_coordinates_byte.nc", time=None, box=False,
                                              coord_type="b", coord_scale=np.float32(0.25),
                                              coords=COORDS / 10)
    out["extras"] = nc_crafted("extras.nc", extras=True)
    out["coords_double"] = nc_crafted("coords_double.nc", coord_type="d", time="d",
                                      times=[0.1, 0.30000001, 0.7, 1.1])
    out["coords_short_f32scale"] = nc_crafted("coords_short_f32scale.nc", coord_type="h",
                                              coord_scale=np.float32(0.01))
    out["coords_int_f64scale"] = nc_crafted("coords_int_f64scale.nc", coord_type="i",
                                            coord_scale=np.float64(0.001), time_scale=np.float64(0.1),
                                            lengths_scale=np.float32(1.1))
    out["coords_float_f64scale"] = nc_crafted("coords_float_f64scale.nc", coord_scale=np.float64(0.3),
                                              time_scale=np.float32(3.3), angles_scale=np.float64(0.7))
    out["scale_near_one"] = nc_crafted("scale_near_one.nc", coord_scale=np.float64(1.0 + 1e-12),
                                       time_scale=np.float32(1.0))
    out["time_irregular"] = nc_crafted("time_irregular.nc", times=[0.1, 0.2, 0.7, 0.71])
    out["time_int"] = nc_crafted("time_int.nc", time="i", times=[3, 5, 9, 10])
    out["time_short_scaled"] = nc_crafted("time_short_scaled.nc", time="h", times=[3, 5, 9, 10],
                                          time_scale=np.float32(0.3))
    out["frame_fixed"] = nc_crafted("frame_fixed.nc", frame_fixed=True)
    out["zero_frames"] = nc_crafted("zero_frames.nc", coords=COORDS[:0])
    out["zero_box"] = nc_crafted("zero_box.nc", lengths_scale=np.float64(0.0))
    # Amber conventions and units
    out["no_conventions"] = nc_crafted("no_conventions.nc", conventions=None)
    out["conventions_other"] = nc_crafted("conventions_other.nc", conventions="CF-1.0")
    out["conventions_list"] = nc_crafted("conventions_list.nc", conventions="CF-1.0,AMBER")
    out["conventions_list_space"] = nc_crafted("conventions_list_space.nc", conventions="CF-1.0, AMBER")
    out["conventions_amber_restart"] = nc_crafted("conventions_amber_restart.nc", conventions="AMBERRESTART")
    out["conventions_lower"] = nc_crafted("conventions_lower.nc", conventions="amber")
    out["no_convention_version"] = nc_crafted("no_convention_version.nc", conv_version=None)
    out["convention_version_2"] = nc_crafted("convention_version_2.nc", conv_version="2.0")
    out["cdf1_bad_conventions"] = nc_crafted("cdf1_bad_conventions.nc", version=1, conventions="X")
    out["spatial_2"] = nc_crafted("spatial_2.nc", spatial=2, coords=COORDS[:, :, :2])
    out["no_spatial"] = nc_crafted("no_spatial.nc", drop=("spatial",), coords=COORDS[:, :, 0])
    out["no_atom"] = nc_crafted("no_atom.nc", drop=("atom",))
    out["no_frame"] = nc_crafted("no_frame.nc", drop=("frame",), time=None, box=False)
    out["no_coordinates"] = nc_crafted("no_coordinates.nc", drop=("coordinates",))
    out["coordinates_no_units"] = nc_crafted("coordinates_no_units.nc", drop=("coordinates.units",))
    out["coordinates_nm"] = nc_crafted("coordinates_nm.nc", units={"coordinates": "nanometer"})
    out["time_ns"] = nc_crafted("time_ns.nc", units={"time": "nanosecond"})
    out["lengths_nm"] = nc_crafted("lengths_nm.nc", units={"cell_lengths": "nm"})
    out["angles_rad"] = nc_crafted("angles_rad.nc", units={"cell_angles": "radian"})
    out["lengths_without_angles"] = nc_crafted("lengths_without_angles.nc", angles=False)
    out["scale_int"] = nc_crafted("scale_int.nc", attrs={"coordinates": {"scale_factor": 2}})
    out["scale_other_variable"] = nc_crafted("scale_other_variable.nc", extras=True,
                                             attrs={"rec_f": {"scale_factor": 2.0}})
    out["scale_array"] = nc_crafted("scale_array.nc",
                                    attrs={"coordinates": {"scale_factor": np.array([1.0, 2.0])}})
    out["no_program"] = nc_crafted("no_program.nc", gattrs={})
    return out


def corrupt_netcdf():
    out = {}
    good = open(tmp("cdf2.nc"), "rb").read()
    out["not_netcdf"] = raw("not_netcdf.nc", b"hello world, not a netcdf file\n")
    out["empty_nc"] = raw("empty_nc.nc", b"")
    out["version_5"] = raw("version_5.nc", b"CDF\x05" + good[4:])
    out["truncated_header"] = raw("truncated_header.nc", good[:100])
    out["truncated_data"] = raw("truncated_data.nc", good[:-40])
    return out


# ---------------------------------------------------------------- Amber ASCII

def mdcrd_text(coords, box=None, title="synthetic mdcrd title", fmt="%8.3f", box_fields=3):
    lines = [title]
    for fr in coords:
        vals = fr.ravel()
        for i in range(0, len(vals), 10):
            lines.append("".join(fmt % v for v in vals[i:i + 10]))
        if box is not None:
            lines.append("".join("%8.3f" % v for v in box[:box_fields]))
    return "\n".join(lines) + "\n"


def rst7_text(coords, title="synthetic restart", time=None, box=None, count=None):
    n = coords.shape[0]
    head = "%5d" % (count if count is not None else n)
    if time is not None:
        head += "%15.7e" % time if isinstance(time, float) else time
    lines = [title, head]
    vals = coords.ravel()
    for i in range(0, len(vals), 6):
        lines.append("".join("%12.7f" % v for v in vals[i:i + 6]))
    if box is not None:
        lines.append("".join("%12.7f" % v for v in box))
    return "\n".join(lines) + "\n"


def amber_files():
    out = {}
    c = COORDS
    t = lambda name, text: raw(name, text.encode("utf-8"))  # noqa: E731
    out["mdcrd_nobox"] = t("nobox.mdcrd", mdcrd_text(c))
    out["mdcrd_box"] = t("box.mdcrd", mdcrd_text(c, box=BOX))
    out["crd_box"] = t("box.crd", mdcrd_text(c, box=BOX))
    out["trj_nobox"] = t("nobox.trj", mdcrd_text(c))
    c10 = (rng.random((3, 10, 3)) * 40 - 20).astype(np.float32)  # 30 values: no short line
    out["mdcrd_10_atoms"] = t("ten.mdcrd", mdcrd_text(c10, box=BOX))
    out["mdcrd_10_atoms_nobox"] = t("ten_nobox.mdcrd", mdcrd_text(c10))
    c4 = (rng.random((3, 4, 3)) * 40 - 20).astype(np.float32)  # 12 values: 2 on the last line
    out["mdcrd_4_atoms_box"] = t("four.mdcrd", mdcrd_text(c4, box=BOX))
    c1 = (rng.random((3, 1, 3)) * 40 - 20).astype(np.float32)
    out["mdcrd_1_atom"] = t("one_atom.mdcrd", mdcrd_text(c1, box=BOX))
    out["mdcrd_crlf"] = raw("crlf.mdcrd", mdcrd_text(c, box=BOX).replace("\n", "\r\n").encode())
    out["mdcrd_trailing_partial"] = t("partial.mdcrd", mdcrd_text(c, box=BOX) + "   1.000   2.000\n")
    out["mdcrd_no_final_newline"] = t("nofinal.mdcrd", mdcrd_text(c).rstrip("\n"))
    out["mdcrd_title_80"] = t("title80.mdcrd", mdcrd_text(c, title="t" * 80 + "   "))
    out["mdcrd_title_81"] = t("title81.mdcrd", mdcrd_text(c, title="t" * 81))
    out["mdcrd_title_unicode"] = t("title_unicode.mdcrd", mdcrd_text(c, title="étude α " * 5))
    out["mdcrd_box_four_fields"] = t("box4.mdcrd", mdcrd_text(c, box=np.array([30, 31, 32, 33.0]),
                                                                box_fields=4))
    out["mdcrd_wide_values"] = t("wide.mdcrd", mdcrd_text(c * 1.0, fmt="%8.2e"))
    out["mdcrd_one_frame_nobox"] = t("one_frame.mdcrd", mdcrd_text(c[:1]))
    out["mdcrd_one_frame_box"] = t("one_frame_box.mdcrd", mdcrd_text(c[:1], box=BOX))
    out["mdcrd_empty"] = t("empty.mdcrd", "")
    out["mdcrd_title_only"] = t("title_only.mdcrd", "title\n")
    bad = mdcrd_text(c).splitlines()
    bad[2] = bad[2][:20] + "  x.yzw " + bad[2][28:]
    out["mdcrd_bad_value"] = t("bad_value.mdcrd", "\n".join(bad) + "\n")
    short = mdcrd_text(c).splitlines()
    short[1] = short[1][:40]
    out["mdcrd_short_line"] = t("short_line.mdcrd", "\n".join(short) + "\n")
    out["mdcrd_stars"] = t("stars.mdcrd", mdcrd_text(c).replace(mdcrd_text(c).splitlines()[1][:8], "********", 1))
    out["rst7_time"] = t("time.rst7", rst7_text(c[0], time=12.5))
    out["rst7_time_box"] = t("time_box.rst7", rst7_text(c[0], time=0.1, box=BOX))
    out["rst7_no_time"] = t("no_time.rst7", rst7_text(c[0]))
    out["rst7_int_time"] = t("int_time.rst7", rst7_text(c[0], time="     7"))
    out["rst7_bad_time"] = t("bad_time.rst7", rst7_text(c[0], time="  abc"))
    out["inpcrd_time"] = t("time.inpcrd", rst7_text(c[0], time=3.0))
    c6 = (rng.random((6, 3)) * 40 - 20).astype(np.float32)
    out["rst7_even"] = t("even.rst7", rst7_text(c6, time=1.0, box=BOX))
    out["rst7_count_mismatch"] = t("count.rst7", rst7_text(c[0], count=N + 1))
    out["rst7_bad_count"] = t("bad_count.rst7", "title\n  abc\n")
    out["rst7_bad_value"] = t("bad_value.rst7", rst7_text(c[0]).replace(rst7_text(c[0]).splitlines()[2][:12],
                                                                       "  bad value ", 1))
    out["rst7_truncated"] = t("truncated.rst7", "\n".join(rst7_text(c[0]).splitlines()[:3]) + "\n")
    return out


# ---------------------------------------------------------------- DCD

def dcd_bytes(coords, *, endian="<", charmm=True, unitcells=None, istart=0, nsavc=1, delta=0.05,
              title=(b"REMARKS crafted dcd".ljust(80),), four_d=False, fixed=None, nset=None,
              marker64=False, extra=None):
    nf, na = coords.shape[:2]

    def rec(payload):
        m = struct.pack(endian + ("q" if marker64 else "i"), len(payload))
        return m + payload + m

    icntrl = [0] * 20
    icntrl[0] = nf if nset is None else nset
    icntrl[1] = istart
    icntrl[2] = nsavc
    icntrl[3] = istart + nsavc * nf
    icntrl[8] = 0 if fixed is None else len(fixed)
    body = bytearray()
    for v in icntrl:
        body += struct.pack(endian + "i", v)
    if charmm:
        body[36:40] = struct.pack(endian + "f", delta)
        icntrl_flags = {10: 1 if unitcells is not None else 0, 11: 1 if four_d else 0, 19: 24}
        for k, v in icntrl_flags.items():
            body[4 * k:4 * k + 4] = struct.pack(endian + "i", v)
    else:
        body[36:44] = struct.pack(endian + "d", delta)
    out = rec(b"CORD" + bytes(body))
    out += rec(struct.pack(endian + "i", len(title)) + b"".join(title))
    out += rec(struct.pack(endian + "i", na))
    free = None
    if fixed is not None:
        free = [i for i in range(na) if i not in fixed]
        out += rec(b"".join(struct.pack(endian + "i", i + 1) for i in free))
    for k, fr in enumerate(coords):
        if extra is not None:
            out += extra
        if charmm and unitcells is not None:
            out += rec(struct.pack(endian + "6d", *unitcells[k % len(unitcells)]))
        sel = fr if (free is None or k == 0) else fr[free]
        for d in range(3):
            out += rec(np.asarray(sel[:, d], dtype=endian + "f4").tobytes())
        if four_d:
            out += rec(np.asarray(sel[:, 0] * 0 + 1, dtype=endian + "f4").tobytes())
    return out


def dcd_files():
    out = {}
    c = COORDS
    out["dcd_mda"] = mda_write("mda.dcd", c, dt=0.5, box=BOX)
    out["dcd_mda_nobox"] = mda_write("mda_nobox.dcd", c, dt=0.5)
    out["dcd_mda_triclinic"] = mda_write("mda_triclinic.dcd", c, box=TRICLINIC, dt=2.0, istart=4, nsavc=2)
    out["dcd_mda_one_frame"] = mda_write("mda_one_frame.dcd", c[:1], box=BOX)
    cos = lambda a: math.cos(math.radians(a))  # noqa: E731
    # [A, gamma, B, beta, alpha, C]
    namd_cos = [[30.5, cos(110.25), 31.25, cos(95.5), cos(80), 32.125]]
    namd_deg = [[30.5, 110.25, 31.25, 95.5, 80.0, 32.125]]
    ortho_cos = [[30.5, 0.0, 31.25, 0.0, 0.0, 32.125], [40.5, 0.0, 41.25, 0.0, 0.0, 42.125]]
    # CHARMM box matrix (lower triangle: a_x, a_y, b_y, a_z, b_z, c_z)
    charmm_matrix = [[30.0, -2.5, 31.0, 1.5, -3.0, 32.0]]
    bad_matrix = [[30.0, 2.5, -31.0, 1.5, 3.0, 32.0]]
    out["dcd_namd_cosines"] = raw("namd_cos.dcd", dcd_bytes(c, unitcells=namd_cos, istart=10, nsavc=5))
    out["dcd_namd_degrees"] = raw("namd_deg.dcd", dcd_bytes(c, unitcells=namd_deg, delta=0.1))
    out["dcd_ortho_cosines"] = raw("ortho.dcd", dcd_bytes(c, unitcells=ortho_cos))
    out["dcd_charmm_matrix"] = raw("charmm_matrix.dcd", dcd_bytes(c, unitcells=charmm_matrix))
    out["dcd_bad_matrix"] = raw("bad_matrix.dcd", dcd_bytes(c, unitcells=bad_matrix))
    out["dcd_big_angles"] = raw("big_angles.dcd", dcd_bytes(c, unitcells=[[30, 190, 31, 90, 90, 32]]))
    out["dcd_zero_cell"] = raw("zero_cell.dcd", dcd_bytes(c, unitcells=[[0, 0, 0, 0, 0, 0]]))
    out["dcd_big_endian"] = raw("big_endian.dcd", dcd_bytes(c, endian=">", unitcells=namd_cos, istart=3,
                                                           nsavc=3, delta=0.002))
    out["dcd_xplor"] = raw("xplor.dcd", dcd_bytes(c, charmm=False, delta=0.0123456789))
    out["dcd_xplor_big_endian"] = raw("xplor_be.dcd", dcd_bytes(c, endian=">", charmm=False, delta=0.5))
    out["dcd_no_cell"] = raw("no_cell.dcd", dcd_bytes(c, istart=7, nsavc=7))
    out["dcd_4d"] = raw("four_d.dcd", dcd_bytes(c, unitcells=ortho_cos, four_d=True))
    out["dcd_fixed"] = raw("fixed.dcd", dcd_bytes(c, unitcells=ortho_cos, fixed={1, 4}))
    out["dcd_wrong_nset"] = raw("wrong_nset.dcd", dcd_bytes(c, nset=99))
    out["dcd_two_titles"] = raw("two_titles.dcd", dcd_bytes(c, title=(b"A".ljust(80), b"B" * 80)))
    out["dcd_no_title"] = raw("no_title.dcd", dcd_bytes(c, title=()))
    out["dcd_marker64"] = raw("marker64.dcd", dcd_bytes(c, unitcells=ortho_cos, marker64=True))
    out["dcd_nsavc_zero"] = raw("nsavc0.dcd", dcd_bytes(c, nsavc=0))
    out["dcd_negative_istart"] = raw("istart_neg.dcd", dcd_bytes(c, istart=-3, nsavc=2))
    full = dcd_bytes(c, unitcells=ortho_cos)
    out["dcd_truncated_frame"] = raw("truncated_frame.dcd", full[:-20])
    out["dcd_extra_bytes"] = raw("extra_bytes.dcd", full + b"\0" * 30)
    header_only = dcd_bytes(c[:0])
    out["dcd_header_only"] = raw("header_only.dcd", header_only)
    out["dcd_truncated_header"] = raw("truncated_header.dcd", full[:50])
    out["dcd_not_cord"] = raw("not_cord.dcd", full[:4] + b"XXXX" + full[8:])
    out["dcd_bad_marker"] = raw("bad_marker.dcd", struct.pack("<i", 85) + full[4:])
    out["dcd_garbage"] = raw("garbage.dcd", b"garbage" * 20)
    out["dcd_empty"] = raw("empty.dcd", b"")
    return out


# ---------------------------------------------------------------- XTC / TRR

def water_like(n_mol, n_frames, seed):
    r = np.random.default_rng(seed)
    out = np.zeros((n_frames, 3 * n_mol, 3), dtype=np.float32)
    for k in range(n_frames):
        o = r.random((n_mol, 3)) * 30
        out[k, 0::3] = o
        out[k, 1::3] = o + [0.9572, 0, 0]
        out[k, 2::3] = o + [-0.24, 0.927, 0]
    return out


def trr_double(coords, box, times):
    nf, na = coords.shape[:2]
    out = b""
    for k in range(nf):
        h = struct.pack(">ii", 1993, 13) + struct.pack(">i", 12) + b"GMX_trn_file"
        sizes = [0, 0, 9 * 8, 0, 0, 0, 0, na * 3 * 8, 0, 0]
        h += struct.pack(">10i", *sizes) + struct.pack(">iii", na, 10 * k, 0)
        h += struct.pack(">dd", times[k], 0.0)
        h += struct.pack(">9d", *box)
        h += np.asarray(coords[k].astype(np.float64) / 10.0, dtype=">f8").tobytes()
        out += h
    return out


def xdr_files():
    out = {}
    c = COORDS
    out["xtc_small"] = mda_write("small.xtc", c, dt=0.5, box=BOX)
    out["trr_small"] = mda_write("small.trr", c, dt=0.5, box=BOX)
    out["xtc_nobox"] = mda_write("nobox.xtc", c, dt=0.1)
    out["trr_nobox"] = mda_write("nobox.trr", c, dt=0.1)
    out["xtc_triclinic"] = mda_write("triclinic.xtc", c, dt=0.3, box=TRICLINIC)
    out["trr_triclinic"] = mda_write("triclinic.trr", c, dt=0.3, box=TRICLINIC)
    out["xtc_one_frame"] = mda_write("one_frame.xtc", c[:1], dt=0.5, box=BOX)
    out["trr_one_frame"] = mda_write("one_frame.trr", c[:1], dt=0.5, box=BOX)
    w = water_like(30, 5, 1)
    out["xtc_water"] = mda_write("water.xtc", w, dt=0.2, box=BOX)
    out["xtc_water_p1000"] = mda_write("water_p1000.xtc", w, dt=0.2, box=TRICLINIC, precision=3)
    out["xtc_water_p2"] = mda_write("water_p2.xtc", w, dt=0.2, box=BOX, precision=2)
    out["xtc_water_p5"] = mda_write("water_p5.xtc", w, dt=0.2, box=BOX, precision=5)
    r = np.random.default_rng(3)
    spread = water_like(14, 3, 5)[:, :40]  # two far clusters: integer sizes beyond 2^24
    spread[:, 21:] += 1.8e6
    out["xtc_large_range"] = mda_write("large_range.xtc", spread, dt=1.0, box=BOX, precision=2)
    wide = water_like(14, 3, 6)[:, :40]  # integer sizes of 24 bits: 72 bits > 64
    wide += (np.repeat(r.random((14, 3)), 3, axis=0)[:40] * 1e5).astype(np.float32)
    out["xtc_wide"] = mda_write("wide.xtc", wide, dt=1.0, box=BOX, precision=3)
    mixed = water_like(20, 3, 4)
    mixed[:, ::7] += 50
    out["xtc_mixed"] = mda_write("mixed.xtc", mixed, dt=1.0, box=BOX)
    out["trr_double"] = raw("double.trr", trr_double(c, [3.05, 0, 0, 0, 3.125, 0, 0, 0, 3.2125],
                                                     [0.0, 0.1, 0.2, 0.3]))
    good_xtc = open(tmp("water.xtc"), "rb").read()
    good_trr = open(tmp("small.trr"), "rb").read()
    out["xtc_truncated"] = raw("truncated.xtc", good_xtc[:-30])
    out["trr_truncated"] = raw("truncated.trr", good_trr[:-30])
    out["xtc_bad_magic"] = raw("bad_magic.xtc", struct.pack(">i", 1996) + good_xtc[4:])
    out["trr_bad_magic"] = raw("bad_magic.trr", struct.pack(">i", 1994) + good_trr[4:])
    out["xtc_empty"] = raw("empty.xtc", b"")
    out["trr_empty"] = raw("empty.trr", b"")
    return out


# ---------------------------------------------------------------- running MDAnalysis

def fhex(x):
    return float(x).hex()


def result_of(u, paths):
    t = u.trajectory
    out = {"n_frames": int(t.n_frames), "dt": fhex(t.dt)}
    frames = []
    for k in range(t.n_frames):
        try:
            ts = t[k]
        except Exception as exc:  # noqa: BLE001
            frames.append({"error": error_of(exc, paths)})
            continue
        dims = ts.dimensions
        pos = np.asarray(ts.positions, dtype=np.float32)
        entry = {"time": fhex(ts.time),
                 "box": None if dims is None else [f32hex(v) for v in np.asarray(dims, dtype=np.float32)]}
        if fhex(t.dt) != out["dt"]:
            entry["dt"] = fhex(t.dt)
        if pos.size > 300:
            entry["pos_sha256"] = hashlib.sha256(pos.astype("<f4").tobytes()).hexdigest()
        else:
            entry["pos"] = pos.astype(">f4").tobytes().hex()
        frames.append(entry)
    out["frames"] = frames
    if type(t).__name__ == "TRJReader":
        out["header"] = t.header
    return out


def error_of(exc, paths):
    message = str(exc)
    for k, p in sorted(enumerate(paths), key=lambda kp: -len(kp[1])):
        message = message.replace(p, "{path%d}" % k)
    return {"type": type(exc).__name__, "message": message}


CASES = []


def case(name, files, fmt=None, dt=None, atoms=N):
    files = [files] if isinstance(files, str) else list(files)
    paths = [tmp(f) for f in files]
    kw = {}
    if fmt is not None:
        kw["format"] = fmt
    if dt is not None:
        kw["dt"] = dt
    entry = {"name": name, "files": files, "format": fmt or "", "dt": None if dt is None else fhex(dt),
             "atoms": atoms}
    try:
        u = mda.Universe.empty(atoms)
        u.load_new(paths[0] if len(paths) == 1 else paths, **kw)
        entry["result"] = result_of(u, paths)
    except Exception as exc:  # noqa: BLE001
        entry["error"] = error_of(exc, paths)
    CASES.append(entry)


def main():
    nc = netcdf_files()
    nc.update(corrupt_netcdf())
    amber = amber_files()
    dcd = dcd_files()
    xdr = xdr_files()

    for name, f in nc.items():
        case(name, f, fmt="NCDF")
        if name in ("mda_box", "mda_one_frame", "no_time", "cdf1", "time_irregular", "coords_double"):
            case(name + "_dt", f, fmt="NCDF", dt=2.5)
    case("mda_box_guessed", nc["mda_box"])
    case("mda_box_atoms_mismatch", nc["mda_box"], fmt="NCDF", atoms=N + 1)
    case("mda_box_format_lower", nc["mda_box"], fmt="ncdf")
    case("mda_box_format_nc", nc["mda_box"], fmt="NC")
    for name, f in amber.items():
        fmt = "INPCRD" if f.endswith((".rst7", ".inpcrd")) else "TRJ"
        atoms = {"mdcrd_10_atoms": 10, "mdcrd_10_atoms_nobox": 10, "mdcrd_4_atoms_box": 4,
                 "mdcrd_1_atom": 1, "rst7_even": 6}.get(name, N)
        case(name, f, fmt=fmt, atoms=atoms)
        if name in ("mdcrd_box", "mdcrd_nobox", "rst7_time", "rst7_no_time"):
            case(name + "_dt", f, fmt=fmt, dt=0.25, atoms=atoms)
    case("mdcrd_guessed", amber["mdcrd_box"])
    case("trj_guessed", amber["trj_nobox"])
    case("crd_guessed", amber["crd_box"])
    case("inpcrd_guessed", amber["inpcrd_time"])
    case("rst7_guessed", amber["rst7_time"])
    case("rst7_atoms_mismatch", amber["rst7_time"], fmt="INPCRD", atoms=N + 1)
    case("mdcrd_atoms_more", amber["mdcrd_box"], fmt="TRJ", atoms=N + 1)
    for name, f in dcd.items():
        case(name, f, fmt="DCD")
        if name in ("dcd_mda", "dcd_namd_cosines"):
            case(name + "_dt", f, fmt="DCD", dt=0.25)
    case("dcd_guessed", dcd["dcd_mda"])
    case("dcd_atoms_mismatch", dcd["dcd_mda"], fmt="DCD", atoms=N + 1)
    for name, f in xdr.items():
        fmt = "XTC" if f.endswith(".xtc") else "TRR"
        atoms = N if name.endswith(("small", "nobox", "triclinic", "one_frame", "double", "trr_truncated",
                                    "trr_bad_magic", "empty")) else None
        if atoms is None:
            atoms = {"xtc_large_range": 40, "xtc_wide": 40, "xtc_mixed": 60}.get(name, 90)
        case(name, f, fmt=fmt, atoms=atoms)
        if name in ("xtc_small", "trr_small", "xtc_one_frame", "trr_one_frame"):
            case(name + "_dt", f, fmt=fmt, dt=0.25, atoms=atoms)
    case("xtc_guessed", xdr["xtc_small"])
    case("trr_guessed", xdr["trr_small"])
    case("xtc_atoms_mismatch", xdr["xtc_small"], fmt="XTC", atoms=N + 1)
    case("trr_atoms_mismatch", xdr["trr_small"], fmt="TRR", atoms=N + 1)

    demo_dir = tempfile.mkdtemp(prefix="demo-", dir=TMP)
    paths = build_demo(demo_dir, n_frames=24, solvent=False)
    shutil.copy(paths["trajectory"], tmp("demo_dry.nc"))

    # chains
    case("chain_nc", [nc["mda_box"], nc["cdf2"], nc["mda_box"]], fmt="NCDF")
    case("chain_nc_dt", [nc["mda_box"], nc["cdf2"]], fmt="NCDF", dt=3.0)
    case("chain_nc_no_time", [nc["no_time"], nc["no_time_no_box"], nc["mda_nobox"]], fmt="NCDF")
    case("chain_nc_no_time_dt", [nc["no_time"], nc["mda_nobox"], nc["no_time"]], fmt="NCDF", dt=0.5)
    case("chain_mdcrd", [amber["mdcrd_box"], amber["mdcrd_nobox"]], fmt="TRJ")
    case("chain_mdcrd_dt", [amber["mdcrd_box"], amber["mdcrd_nobox"]], fmt="TRJ", dt=0.75)
    case("chain_rst7", [amber["rst7_time"], amber["rst7_no_time"], amber["rst7_time_box"]], fmt="INPCRD")
    case("chain_rst7_dt", [amber["rst7_no_time"], amber["rst7_time"], amber["rst7_no_time"]], fmt="INPCRD",
         dt=4.0)
    case("chain_dcd", [dcd["dcd_mda"], dcd["dcd_namd_cosines"]], fmt="DCD")
    case("chain_dcd_dt", [dcd["dcd_mda"], dcd["dcd_namd_cosines"]], fmt="DCD", dt=0.5)
    case("chain_xdr_guessed", [xdr["xtc_small"], xdr["trr_small"], xdr["trr_triclinic"]])
    case("chain_xdr_guessed_dt", [xdr["xtc_small"], xdr["trr_small"]], dt=3.0)
    case("chain_mixed", [nc["mda_box"], nc["no_time"], dcd["dcd_mda"], xdr["xtc_small"],
                         amber["inpcrd_time"], amber["trj_nobox"], amber["mdcrd_box"]])
    case("chain_mixed_dt", [amber["trj_nobox"], nc["no_time"], dcd["dcd_namd_cosines"], nc["mda_box"]],
         dt=0.4)
    case("chain_rst7_guessed", [nc["mda_box"], amber["rst7_time"]])
    case("chain_crd_guessed", [nc["mda_box"], amber["crd_box"]])
    case("chain_bad_second", [nc["mda_box"], nc["no_conventions"]], fmt="NCDF")
    case("chain_bad_first", [nc["cdf1"], nc["mda_box"]], fmt="NCDF")
    case("chain_value_error_first", [amber["rst7_bad_time"], nc["mda_box"]])
    case("chain_dcd_error", [nc["mda_box"], dcd["dcd_garbage"]])
    case("chain_atoms_mismatch", [xdr["xtc_small"], xdr["trr_small"]], atoms=N + 1)
    case("chain_atoms_differ", [xdr["xtc_small"], xdr["xtc_water"]], atoms=N)
    case("chain_atoms_differ_first", [xdr["xtc_water"], xdr["xtc_small"], xdr["trr_small"]], atoms=N)
    case("chain_atoms_differ_nc", [nc["mda_box"], "demo_dry.nc"], atoms=N)
    case("chain_truncated", [xdr["xtc_truncated"], xdr["xtc_water"]], atoms=90)
    case("chain_single", [nc["mda_box"]], fmt="NCDF")

    # the demo trajectory (MDATK's demo_dry fixture)
    case("demo_dry", "demo_dry.nc", fmt="NCDF", atoms=773)
    case("demo_dry_dt", "demo_dry.nc", fmt="NCDF", dt=5.0, atoms=773)
    demo = mda.Universe(paths["topology"], paths["trajectory"], format="NCDF")
    demo_xyz = np.array([ts.positions.copy() for ts in demo.trajectory])
    demo_box = demo.trajectory.ts.dimensions
    mda_write("demo_dry.xtc", demo_xyz[:6], dt=10.0, box=demo_box)
    mda_write("demo_dry.trr", demo_xyz[:3], dt=10.0, box=demo_box)
    mda_write("demo_dry.dcd", demo_xyz[:3], dt=10.0, box=demo_box)
    case("demo_dry_xtc", "demo_dry.xtc", fmt="XTC", atoms=773)
    case("demo_dry_trr", "demo_dry.trr", fmt="TRR", atoms=773)
    case("demo_dry_dcd", "demo_dry.dcd", fmt="DCD", atoms=773)

    # the NetCDF writer: MDAnalysis' output for the frames of these cases (the test writes the
    # frames it reads from the source file and compares bytes)
    writer_cases = []
    for name, src, n_frames in [("demo_dry", "demo_dry.nc", None), ("box", nc["mda_triclinic"], None),
                                ("nobox", nc["mda_nobox"], None), ("one_frame", nc["mda_one_frame"], None),
                                ("from_xtc", "demo_dry.xtc", 2), ("dcd", dcd["dcd_namd_cosines"], None)]:
        u = mda.Universe.empty(773 if src.startswith("demo") else N)
        u.load_new(tmp(src))
        frames = list(range(u.trajectory.n_frames if n_frames is None else n_frames))
        out = os.path.join(WRITER, name + ".nc")
        with mda.Writer(out, n_atoms=u.atoms.n_atoms, format="NCDF") as w:
            for k in frames:
                u.trajectory[k]
                w.write(u.atoms)
        expected = "writer/" + name + ".nc"
        if open(out, "rb").read() == open(tmp(src), "rb").read():
            os.remove(out)  # the source file itself (MDATK wrote it with NCDFWriter)
            expected = src
        writer_cases.append({"name": name, "source": src, "frames": frames, "atoms": int(u.atoms.n_atoms),
                             "expected": expected})

    # copy the data files (not MDAnalysis' offset caches)
    used = set()
    for c in CASES:
        used.update(c["files"])
    for f in sorted(used):
        if f == "demo_dry.nc":
            shutil.copy(tmp(f), os.path.join(OUT, f))
        else:
            shutil.copy(tmp(f), os.path.join(FILES, f))
    for c in CASES:
        c["files"] = [f if f == "demo_dry.nc" else "files/" + f for f in c["files"]]
    for w in writer_cases:
        for key in ("source", "expected"):
            if not w[key].startswith(("demo_dry.nc", "writer/")):
                w[key] = "files/" + w[key]

    with open(os.path.join(OUT, "cases.json"), "w") as f:  # one case per line
        f.write('{"cases": [\n')
        f.write(",\n".join(json.dumps(c, ensure_ascii=False) for c in CASES))
        f.write('\n],\n"writer": [\n')
        f.write(",\n".join(json.dumps(w, ensure_ascii=False) for w in writer_cases))
        f.write("\n]}\n")

    formats = {"mdatk_table": loader._TRAJECTORY_FORMATS, "guess": []}
    for name in ["a.nc", "a.NC", "a.ncdf", "a.netcdf", "a.mdcrd", "a.crd", "a.trj", "a.rst7", "a.inpcrd",
                 "a.restrt", "a.ncrst", "a.dcd", "a.DCD", "a.xtc", "a.trr", "a.pdb", "a.ent", "a.gro",
                 "a.xyz", "noext", "dir.d/noext", "a.b.Xtc", ".xtc", "a.", "a.tng", "a.mdcrd.gz", "a.nc.bz2"]:
        path = tmp(name)
        os.makedirs(os.path.dirname(path), exist_ok=True)
        open(path, "w").close()
        entry = {"path": name}
        try:
            entry["guess"] = mda.lib.util.guess_format(name)
        except Exception as exc:  # noqa: BLE001
            entry["guess_error"] = error_of(exc, [name])
        try:
            u = mda.Universe.empty(N)
            u.load_new(path)
        except Exception as exc:  # noqa: BLE001
            entry["load_error"] = error_of(exc, [path])
        formats["guess"].append(entry)
    for fmt in ["NCRST", "RST7", "NETCDF", "CRD", "XYZ", "bogus"]:
        try:
            u = mda.Universe.empty(N)
            u.load_new(tmp(nc["mda_box"]), format=fmt)
            formats.setdefault("format_errors", []).append({"format": fmt, "loaded": True})
        except Exception as exc:  # noqa: BLE001
            formats.setdefault("format_errors", []).append(
                {"format": fmt, "error": error_of(exc, [tmp(nc["mda_box"])])})
    with open(os.path.join(OUT, "formats.json"), "w") as f:
        json.dump(formats, f, indent=1, ensure_ascii=False)
        f.write("\n")
    shutil.rmtree(TMP, ignore_errors=True)


main()
