"""MDBridge's improvements of the Amber / DCD / XTC / TRR readers, applied to the expected data.

MDBridge fixes the reader bugs of MDAnalysis 2.11.0 listed in docs/IMPROVEMENTS.md, sections 2.1
(parm7), 2.4 (AMBER NetCDF), 2.5 (AMBER ASCII mdcrd / rst7) and 2.6 (DCD / XTC / TRR); each fix is
described in docs/DIFFERENCES.md ("Readers and topologies"). The generators still record what
MDAnalysis does; this module turns that into what MDBridge does, on the test side (pure Python and
numpy, plus MDAnalysis' public API where MDAnalysis can compute the improved result from a repaired
copy of a file). Never MDAnalysis code: only behaviour.

    import improvements_readers_amber
    improvements_readers_amber.apply()     # independent_tables.apply() calls it

* apply(): RD-P12 (and MDBridge's analysis versions, energy 2 among them, given to MDATK's analyses
  by analysis_versions.apply_mdatk() when MDATK is importable). MDAnalysis' TOPParser rounds each CHARGE field to float32 and divides by
  float32(18.2223) in float32; MDBridge divides the field by 18.2223 in float64. The TOPParser's
  parse() is wrapped (public API: the returned Topology's `charges` values are replaced) so that
  every Universe built from a parm7 file - topology dumps, selections, analyses (energy), the
  command line - carries MDBridge's charges.
* apply() also gives MDATK's parm7 section reader (mdatk.io.parm7.read_sections, which the energy
  analysis' LJ table and make_parm7_cases.py use) MDATK 0.1.0's format lines (same owner; its code
  ported): '(a80)' without a repeat count (CHAMBER's CTITLE), '(8(F9.5))' (ff19SB's CMAP
  parameters), and a section with a format it cannot parse is skipped when not wanted (CHAMBER's
  '(i2,a78)'). MDATK 0.0.4.1 raised "cannot parse format line" for all three.
* apply() also gives MDATK version 0.1.0's time-step override (LD-05; io/loader.py _FrameTimes,
  pipeline.py input_fingerprint): with inputs.dt_ps set, load_universe's trajectory has frame k at
  k x dt_ps whatever the files' time stamps, and the cache keys' input fingerprint gets
  "time_axis": "frame x dt_ps". adjust_trajectory_cases() gives the trajectory cases with a dt the
  same times (MDBridge's openTrajectories applies the override).
* adjust_parm7_topologies(...): make_parm7_cases.py's topology variants (RD-P01 ... RD-P14):
  the error texts MDBridge gives, and the topologies of the files MDBridge now reads (dumped by
  MDAnalysis from a repaired copy).
* adjust_trajectory_cases(...): make_trajectory_cases.py's cases (RD-N*, RD-A*, RD-D*), with the
  reader warnings MDBridge issues for each ("warnings", LD-03; trajectory_warnings()).
* netcdf_extra_files(...): the crafted NetCDF files of RD-N05 (streaming record count) and RD-N07
  (cell variables without the frame dimension), which MDAnalysis cannot read correctly.
"""

from __future__ import annotations

import json
import math
import os
import re

import numpy as np

AMBER_CHARGE_FACTOR = 18.2223

# ---------------------------------------------------------------------------------------- RD-P12

_FORMAT_RE = re.compile(r"\(\s*(\d*)\s*([aAiIeEfF])\s*(\d+)")


def parm7_charges(path):
    """The CHARGE fields of a parm7 file divided by 18.2223 in float64 (None without CHARGE)."""
    with open(path, encoding="utf-8", errors="replace") as fh:
        lines = fh.read().splitlines()
    natom = None
    for k, line in enumerate(lines):
        if line.startswith("%FLAG POINTERS"):
            j = k + 1
            while lines[j].startswith(("%FORMAT", "%COMMENT")):
                j += 1
            natom = int(lines[j].split()[0])
            break
    for k, line in enumerate(lines):
        if not line.startswith("%FLAG CHARGE"):
            continue
        j = k + 1
        while lines[j].startswith("%COMMENT"):
            j += 1
        m = _FORMAT_RE.search(lines[j])
        per, width = int(m.group(1) or 1), int(m.group(3))
        values = []
        for data in lines[j + 1:j + 1 + math.ceil(natom / per)]:
            for c in range(0, min(len(data), per * width), width):
                text = data[c:c + width].strip()
                if text:
                    values.append(float(text.replace("D", "E").replace("d", "e")) / AMBER_CHARGE_FACTOR)
        return np.array(values, dtype=np.float64)
    return None


def apply():
    """Wrap TOPParser.parse so that charges are MDBridge's float64 ones; MDATK's analyses get
    MDBridge's versions (idempotent)."""
    from MDAnalysis.topology.TOPParser import TOPParser

    # the energy results change with the charges: MDBridge's analysis versions (energy 2), from
    # the one table of them (analysis_versions.py), so the order of the improvement modules is moot
    _sibling("analysis_versions").apply_mdatk()
    _patch_mdatk_parm7()
    _patch_mdatk_time_step()
    if getattr(TOPParser, "_mdbridge_float64_charges", False):
        return
    original = TOPParser.parse

    def parse(self, **kwargs):
        top = original(self, **kwargs)
        charges = parm7_charges(self.filename)
        if charges is not None and hasattr(top, "charges") and len(charges) == len(top.charges.values):
            top.charges.values = charges
        return top

    TOPParser.parse = parse
    TOPParser._mdbridge_float64_charges = True


# MDATK 0.1.0, mdatk/io/parm7.py (_FMT and read_sections)
_FMT_010 = re.compile(r"\(\s*(\d*)\s*\(?\s*([aAiIeEfF])\s*(\d+)(?:\.(\d+))?\s*\)?\s*\)")


def _patch_mdatk_parm7():
    """MDATK's parm7 section reader with MDATK 0.1.0's format lines (idempotent; a no-op without
    MDATK)."""
    try:
        from mdatk.io import parm7
    except ImportError:
        return
    if getattr(parm7, "_mdbridge_format_lines", False):
        return
    Section = parm7.Section

    def read_sections(path, wanted=None):
        wanted_set = set(wanted) if wanted is not None else None
        out = {}
        current = None

        def flush():
            if current is not None and (wanted_set is None or current.name in wanted_set):
                out[current.name] = current.values()

        with open(path, "r", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if line.startswith("%FLAG"):
                    flush()
                    current = Section(line[5:].strip(), 0, "a", 80, [])
                    continue
                if line.startswith("%FORMAT"):
                    if current is None:
                        continue
                    m = _FMT_010.search(line)
                    if not m:
                        if wanted_set is not None and current.name not in wanted_set:
                            current = None  # a section nobody asked for (CHAMBER's '(i2,a78)'): not read
                            continue
                        raise ValueError(f"{path}: cannot parse format line {line.strip()!r}")
                    current.count = int(m.group(1) or 1)
                    current.kind = m.group(2).lower()
                    current.width = int(m.group(3))
                    continue
                if line.startswith("%VERSION") or line.startswith("%COMMENT"):
                    continue
                if current is not None and (wanted_set is None or current.name in wanted_set):
                    current.lines.append(line)
        flush()
        return out

    parm7._FMT = _FMT_010
    parm7.read_sections = read_sections
    parm7._mdbridge_format_lines = True


class _FrameTimes:
    """MDATK 0.1.0 (io/loader.py): trajectory transformation of the time-step override, the time of
    frame k is k * dt."""

    def __init__(self, dt_ps):
        self.dt_ps = dt_ps

    def __call__(self, ts):
        ts.time = ts.frame * self.dt_ps
        return ts


def _patch_mdatk_time_step():
    """LD-05: MDATK's load_universe with MDATK 0.1.0's override of the files' times, and its input
    fingerprint with the "time_axis" entry (idempotent; a no-op without MDATK)."""
    try:
        from mdatk.io import loader
    except ImportError:
        return
    if getattr(loader, "_mdbridge_frame_times", False):
        return
    original = loader.load_universe

    def load_universe(topology, trajectories=(), dt_ps=None):
        u = original(topology, trajectories, dt_ps)
        if dt_ps and dt_ps == dt_ps:  # MDBridge: also coordinates read from the topology file
            u.trajectory.add_transformations(_FrameTimes(float(dt_ps)))
        return u

    loader.load_universe = load_universe
    loader._mdbridge_frame_times = True
    import sys

    pipeline = sys.modules.get("mdatk.pipeline")
    if pipeline is None:
        from mdatk import pipeline
    pipeline.load_universe = load_universe
    fingerprint = pipeline.Pipeline.input_fingerprint

    def input_fingerprint(self, inputs):
        fp = fingerprint(self, inputs)
        if self.project.inputs.dt_ps:  # results cached with an override before it set the time axis
            fp["time_axis"] = "frame x dt_ps"
        return fp

    pipeline.Pipeline.input_fingerprint = input_fingerprint


def _sibling(name):
    """A module of tests/parity (this folder), whether or not it is on sys.path."""
    import importlib.util
    import sys
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                                                     name + ".py"))
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


# ------------------------------------------------------------------------- parm7 topology variants

def _line(text, flag, offset):
    """1-based number of the line `offset` lines after the line `%FLAG <flag>` of `text`."""
    lines = re.split(r"\r\n|\r|\n", text)
    return lines.index("%FLAG " + flag) + 1 + offset


def _replace_section_line(text, flag, data_index, new):
    lines = text.split("\n")
    k = lines.index("%FLAG " + flag) + 2 + data_index
    lines[k] = new
    return "\n".join(lines)


def adjust_parm7_topologies(folder, variants, errors, small_text, dump):
    """Rewrites make_parm7_cases.py's topology outputs to MDBridge's behaviour.

    folder: where the variants and their dumps are; variants: [(name, text)]; errors: the list of
    {"file", "stage", "error"} MDAnalysis gave; small_text: the unmodified small system; dump: the
    dump function (dump_topologies.dump). Returns the new error list (in variant order)."""
    texts = dict(variants)
    path = "{path}"
    value_error = lambda message: {"type": "ValueError", "message": message}  # noqa: E731
    mass_lines = math.ceil(40 / 5)
    new_errors = {
        # RD-P05: a blank field before the last non-blank one
        "mass_blank_field": value_error(f"{path} line {_line(texts['mass_blank_field'], 'MASS', 2)}: "
                                        "blank field 1 in %FLAG MASS"),
        # RD-P03: a second line of '%FLAG' alone is the header check's own message
        "title_one_token": value_error(f"{path} is not a valid TOP file. 'TITLE' missing in header"),
        # RD-P02 / RD-P04
        "no_pointers": value_error(f"{path}: required section %FLAG POINTERS is missing"),
        "not_flag_after_section": value_error(
            f"{path} line {_line(texts['not_flag_after_section'], 'MASS', 2 + mass_lines)}: expected a %FLAG "
            f"line after the {mass_lines} data lines of %FLAG MASS, found 'trailing junk'"),
        "empty_line_after_section": value_error(
            f"{path} line {_line(texts['empty_line_after_section'], 'MASS', 2 + mass_lines)}: expected a %FLAG "
            f"line after the {mass_lines} data lines of %FLAG MASS, found ''"),
        # RD-P09
        "atomic_number_119": value_error(f"{path}: ATOMIC_NUMBER 119 of atom 1 is not an element"),
        # RD-P04
        "no_atom_name": value_error(f"{path}: required section %FLAG ATOM_NAME is missing"),
        "no_angles": value_error(f"{path}: required section %FLAG ANGLES_WITHOUT_HYDROGEN is missing"),
        # RD-P02 (RM): a misspelt %FORMAT line
        "bad_format_line": value_error(
            f"{path} line {_line(texts['bad_format_line'], 'MASS', 1)}: expected %FORMAT or %COMMENT after "
            "%FLAG MASS, found '%FROMAT(5E16.8)'"),
        # RD-P03 (RM): an integer field
        "bad_int": value_error(f"{path} line {_line(texts['bad_int'], 'RESIDUE_POINTER', 2)}, "
                               "%FLAG RESIDUE_POINTER field 2: '1x' is not an integer"),
        # RD-P08: RESIDUE_POINTER validated
        "resptr_zero": value_error(f"{path}: %FLAG RESIDUE_POINTER residue 1 starts at atom 0; the first "
                                   "residue must start at atom 1"),
        "resptr_unsorted": value_error(f"{path}: %FLAG RESIDUE_POINTER is not increasing at residue 3 "
                                       "(11, then 6)"),
        "residue_pointer_offset": value_error(f"{path}: %FLAG RESIDUE_POINTER residue 1 starts at atom 3; "
                                              "the first residue must start at atom 1"),
    }
    # files MDBridge reads, as MDAnalysis reads a repaired copy
    repaired = {
        "unknown_first": small_text,          # RD-P01: the unknown section is skipped
        "pointer_lines_short": small_text,    # RD-P03: 30 POINTERS values are enough
        "charge_three_per_line": small_text,  # RD-P06: lines counted from %FORMAT(3E24.16)
        # RD-P05: a blank chain ID before others is the chain "" (MDAnalysis drops it); repaired
        # with the chain "Z", which is then renamed ""
        "chain_blank": texts["chain_blank"].replace("A   A       B   ", "A   A   Z   B   "),
    }
    for name, text in repaired.items():
        tmp = os.path.join(folder, name + ".repaired.parm7")
        with open(tmp, "w", encoding="utf-8", newline="") as fh:
            fh.write(text)
        data = dump(tmp)
        os.remove(tmp)
        data["file"] = name + ".parm7"
        if name == "chain_blank":
            if "Z" not in data["atoms"]["chainIDs"]:
                raise RuntimeError("chain_blank: the repaired chain is missing")
            data["atoms"]["chainIDs"] = ["" if c == "Z" else c for c in data["atoms"]["chainIDs"]]
            data["segments"]["segids"] = ["" if s == "Z" else s for s in data["segments"]["segids"]]
        with open(os.path.join(folder, name + ".parm7.topology.json"), "w", encoding="utf-8") as fh:
            json.dump(data, fh, indent=1, ensure_ascii=False)
            fh.write("\n")
    for name in new_errors:
        dumped = os.path.join(folder, name + ".parm7.topology.json")
        if os.path.exists(dumped):
            os.remove(dumped)
    old = {e["file"]: e for e in errors}
    out = []
    for name, _text in variants:
        if name in repaired:
            continue
        if name in new_errors:
            out.append({"file": name + ".parm7", "stage": "parser", "error": new_errors[name]})
        elif name + ".parm7" in old:
            out.append(old[name + ".parm7"])
    return out


# ------------------------------------------------------------------------------- NetCDF helpers

def convert_netcdf(src, dst, *, drop=(), cell_per_frame=False):
    """A CDF-2 copy of a NetCDF file (scipy), without the variables `drop`; with cell_per_frame,
    cell variables without the frame dimension become (frame, ...) copies."""
    import scipy.io

    with scipy.io.netcdf_file(src, "r", mmap=False) as f:
        out = scipy.io.netcdf_file(dst, "w", version=2)
        for k, v in f._attributes.items():
            setattr(out, k, v)
        nframes = f.variables["coordinates"].shape[0]
        for name, size in f.dimensions.items():
            out.createDimension(name, size)
        for name, var in f.variables.items():
            if name in drop:
                continue
            dims = var.dimensions
            data = var[:].copy()
            if cell_per_frame and name.startswith("cell_") and "frame" not in dims and len(dims) == 1 \
                    and dims[0] in ("cell_spatial", "cell_angular"):
                dims = ("frame",) + dims
                data = np.repeat(data[None, :], nframes, axis=0)
            new = out.createVariable(name, var.typecode(), dims)
            for k, v in var._attributes.items():
                setattr(new, k, v)
            if data.shape:
                new[:] = data
            else:
                new.assignValue(data)
        out.close()


def _record_layout(path):
    """(record begin, record stride, records) of a NetCDF classic file, from its header (scipy)."""
    import scipy.io

    with scipy.io.netcdf_file(path, "r", mmap=False) as f:
        records = f.variables["coordinates"].shape[0]
        sizes = []
        for var in f.variables.values():
            if var.isrec:
                n = var.itemsize() * int(np.prod(var.shape[1:]))
                sizes.append(n + ((4 - n % 4) % 4 if var.typecode() in "bch" else 0))
        stride = sizes[0] if len(sizes) == 1 else sum(sizes)
    return os.path.getsize(path) - records * stride, stride, records


def netcdf_extra_files(tmp, raw, nc_crafted, coords):
    """The crafted files of RD-N05 / RD-N07: names -> file names (in the generator's TMP)."""
    import scipy.io

    out = {}
    good = open(tmp("cdf2.nc"), "rb").read()
    out["streaming"] = raw("streaming.nc", good[:4] + b"\xff\xff\xff\xff" + good[8:])
    for name, dims in (("cell_constant", ("cell_spatial",)), ("cell_per_frame_scalar", ("frame",))):
        nc_crafted(name + ".nc", box=False, coords=coords)
        f = scipy.io.netcdf_file(tmp(name + ".nc"), "a", version=2)
        lv = f.createVariable("cell_lengths", "d", dims)
        lv.units = "angstrom"
        lv[:] = np.array([30.1, 31.7, 32.3] if dims == ("cell_spatial",) else [30.1 + k for k in range(len(coords))])
        av = f.createVariable("cell_angles", "d", ("cell_angular",))
        av.units = "degree"
        av[:] = np.array([90.0, 91.3, 119.9])
        f.close()
        out[name] = name + ".nc"
    return out


# --------------------------------------------------------------------------- trajectory cases

def _f32bits(values):
    return [np.asarray(v, dtype=">f4").tobytes().hex() for v in values]


def adjust_trajectory_cases(cases, tmp, rerun):
    """Rewrites make_trajectory_cases.py's CASES (in place) to MDBridge's behaviour.

    tmp(name): the path of a generated file; rerun(files, fmt, dt, atoms, paths_for_messages):
    MDAnalysis' result entry ({"result": ...} or {"error": ...}) for files given as paths."""
    by_name = {c["name"]: c for c in cases}

    def set_error(name, type_, message):
        c = by_name[name]
        c.pop("result", None)
        c["error"] = {"type": type_, "message": message}

    def set_result(name, result):
        c = by_name[name]
        c.pop("error", None)
        c["result"] = result

    def path_of(name):
        return tmp(by_name[name]["files"][0])

    def dt_of(c):
        return None if c["dt"] is None else float.fromhex(c["dt"])

    # RD-N08: CDF-1 files are read like CDF-2 ones
    for name in ("cdf1", "cdf1_dt", "chain_bad_first"):
        c = by_name[name]
        paths = []
        for k, f in enumerate(c["files"]):
            src = tmp(f)
            if open(src, "rb").read(4) == b"CDF\x01":
                dst = tmp(f"cdf1_as_cdf2_{k}.nc")
                convert_netcdf(src, dst)
                src = dst
            paths.append(src)
        entry = rerun(paths, c["format"] or None, dt_of(c), c["atoms"])
        set_result(name, entry["result"])
    # RD-N01: lengths (0, 0, 0) are no cell (MDAnalysis reads the file without its cell)
    dst = tmp("zero_box_without_cell.nc")
    convert_netcdf(path_of("zero_box"), dst, drop=("cell_lengths", "cell_angles"))
    set_result("zero_box", rerun([dst], "NCDF", None, by_name["zero_box"]["atoms"])["result"])
    # RD-N02: the atom-count message of the other readers
    set_error("mda_box_atoms_mismatch", "ValueError",
              "The topology and NCDF trajectory files don't have the same number of atoms!\n"
              "Topology number of atoms 8\nTrajectory: {path0} Number of atoms 7")
    set_error("chain_atoms_differ_nc", "ValueError",
              "The following trajectories do not have the correct n_atoms  (7):\n['{path1}']")
    # RD-N03
    set_error("coordinates_no_units", "OSError",
              "{path0}: variable 'coordinates' has no 'units' attribute (expected 'angstrom')")
    # RD-N04
    set_error("empty_nc", "OSError", "{path0}: file is empty")
    set_error("version_5", "OSError", "{path0}: unsupported NetCDF variant (CDF-5)")
    set_error("truncated_header", "OSError", "{path0}: truncated header")
    begin, stride, records = _record_layout(path_of("cdf2"))
    size = os.path.getsize(path_of("truncated_data"))
    set_error("truncated_data", "OSError",
              f"{{path0}}: file ends inside frame {(size - begin) // stride} (expected "
              f"{begin + records * stride} bytes, found {size})")
    # RD-N05
    set_error("zero_frames", "OSError", "{path0}: trajectory has no frames")
    set_result("streaming", json.loads(json.dumps(by_name["cdf2"]["result"])))
    # RD-N06: dt from the scaled times (as the frames' times)
    import scipy.io
    for name in ("mda_scaled", "coords_int_f64scale", "coords_float_f64scale", "time_short_scaled",
                 "scale_near_one"):
        c = by_name[name]
        with scipy.io.netcdf_file(path_of(name), "r", mmap=False) as f:
            t = f.variables["time"]
            raw = np.asarray(t[:2]).copy()
            factor = t._attributes.get("scale_factor")
        if factor is not None:
            scale = float(np.asarray(factor).ravel()[0])
            if not abs(scale - 1.0) <= 1e-9 * max(abs(scale), 1.0):
                raw = raw * np.asarray(factor).ravel()[0]  # numpy's promotion, as the frames' times
        c["result"]["dt"] = float(raw[1] - raw[0]).hex()
    # RD-N07: a cell variable without the frame dimension is every frame's cell
    dst = tmp("cell_constant_per_frame.nc")
    convert_netcdf(path_of("cell_constant"), dst, cell_per_frame=True)
    set_result("cell_constant", rerun([dst], "NCDF", None, by_name["cell_constant"]["atoms"])["result"])
    set_error("cell_per_frame_scalar", "OSError",
              "{path0}: variable 'cell_lengths' has dimensions (frame); expected (frame, cell_spatial) "
              "or (cell_spatial)")

    # RD-A01: one frame without a box line (MDAnalysis reads a copy with the frame twice)
    src = path_of("mdcrd_one_frame_nobox")
    text = open(src, encoding="utf-8").read().splitlines()
    dst = tmp("one_frame_twice.mdcrd")
    with open(dst, "w", encoding="utf-8") as fh:
        fh.write("\n".join(text + text[1:]) + "\n")
    result = rerun([dst], "TRJ", None, by_name["mdcrd_one_frame_nobox"]["atoms"])["result"]
    result["n_frames"], result["frames"] = 1, result["frames"][:1]
    set_result("mdcrd_one_frame_nobox", result)
    # RD-A03: one atom with box lines (MDAnalysis reads each box line as a frame)
    c = by_name["mdcrd_1_atom"]
    frames = c["result"]["frames"]
    out = []
    for k in range(0, len(frames), 2):
        box = np.frombuffer(bytes.fromhex(frames[k + 1]["pos"]), dtype=">f4")
        out.append({"time": float(k // 2).hex(),
                    "box": _f32bits(list(box) + [90.0, 90.0, 90.0]), "pos": frames[k]["pos"]})
    c["result"]["n_frames"], c["result"]["frames"] = len(out), out
    # RD-A02 / RD-A03: the layout is checked against the atom count
    set_error("mdcrd_box_four_fields", "OSError",
              "{path0} line 5: expected a box line (3 values) or the first 10 coordinates of the next "
              "frame, found 4 values (does the topology match?)")
    wide = open(path_of("mdcrd_wide_values"), encoding="utf-8").read().splitlines()[1].rstrip()
    set_error("mdcrd_wide_values", "OSError",
              f"{{path0}} line 2: expected 10 coordinates for 7 atoms, found {math.ceil(len(wide) / 8)} "
              "(does the topology match?)")
    set_error("mdcrd_short_line", "OSError",
              "{path0} line 2: expected 10 coordinates for 7 atoms, found 5 (does the topology match?)")
    set_error("mdcrd_atoms_more", "OSError",
              "{path0} line 4: expected 4 coordinates for 8 atoms, found 1 (does the topology match?)")
    # RD-A05
    bad = open(path_of("mdcrd_bad_value"), encoding="utf-8").read().splitlines()[2]
    for k in range(10):  # the first 8-column field that is not a number
        field = bad[8 * k:8 * k + 8]
        try:
            float(field)
        except ValueError:
            break
    set_error("mdcrd_bad_value", "OSError",
              f"{{path0}} line 3, columns {8 * k + 1}-{8 * k + 8}: {field!r} is not a number")
    set_error("mdcrd_stars", "OSError", "{path0} line 2, columns 1-8: '********' is not a number "
                                        "(Fortran overflow: the value does not fit the field)")
    set_error("rst7_bad_time", "OSError", "{path0} line 2: time 'abc' is not a number")
    set_error("rst7_count_mismatch", "OSError",
              "{path0} line 6: expected 6 coordinates, found 3 (the file has 8 atoms)")
    set_error("rst7_bad_count", "OSError", "{path0} line 2: atom count 'abc' is not a non-negative integer")
    set_error("rst7_bad_value", "OSError", "{path0} line 3, columns 1-12: '  bad value ' is not a number")
    set_error("rst7_truncated", "OSError", "{path0}: the file ends at line 3; 7 atoms need 6 lines")
    # RD-A06
    set_error("mdcrd_empty", "OSError", "{path0}: file is empty")
    set_error("mdcrd_title_only", "OSError", "{path0}: no frames after the title")

    # RD-D01: NSAVC 0 is read as 1 (MDAnalysis reads a copy with NSAVC 1)
    data = bytearray(open(path_of("dcd_nsavc_zero"), "rb").read())
    data[16:20] = (1).to_bytes(4, "little")
    dst = tmp("nsavc1.dcd")
    open(dst, "wb").write(bytes(data))
    set_result("dcd_nsavc_zero", rerun([dst], "DCD", None, by_name["dcd_nsavc_zero"]["atoms"])["result"])
    # RD-D02
    set_error("dcd_4d", "OSError", "{path0}: 4-dimensional DCD files are not supported")
    # RD-D03: XDR errors name the file and the frame
    for name in ("xtc_truncated", "trr_truncated", "chain_truncated"):
        for k, frame in enumerate(by_name[name]["result"]["frames"]):
            if "error" in frame:  # the first file's frame k in each of these cases
                frame["error"]["message"] = f"{{path0}} frame {k}: " + frame["error"]["message"]
    for name in ("xtc_bad_magic", "xtc_empty", "trr_empty"):
        by_name[name]["error"]["message"] = "{path0}: " + by_name[name]["error"]["message"]
    set_error("rst7_empty", "OSError", "{path0}: file is empty")  # RD-B06

    # LD-03: the warnings of the cases MDBridge reads (RD-N08, RD-N10, RD-A04, RD-D01, LD-05)
    single = {(c["files"][0], c["dt"]): c["result"] for c in cases if len(c["files"]) == 1 and "result" in c}
    for c in cases:
        if "result" in c:
            found = trajectory_warnings(c, tmp, single)
            if found:
                c["warnings"] = found

    # LD-05 (MDATK 0.1.0): a given dt replaces the files' times, frame k at k x dt (after the warnings
    # above, which compare the files' own steps)
    for c in cases:
        dt = float.fromhex(c["dt"]) if c.get("dt") is not None else 0.0
        if "result" in c and dt:
            for k, frame in enumerate(c["result"]["frames"]):
                if "time" in frame:
                    frame["time"] = (k * dt).hex()


NO_TIMES = ": the file has no time information; frames are taken to be 1 ps apart (set inputs.dt_ps for the " \
           "real time step)"


def _time_step_warning(c, k, f, stored, single):
    """LD-05 (I3, src/io/trajectory.cpp): inputs.dt_ps given and file k's own step (frame 1 - frame
    0 as read) differs from it. single: (file, dt) -> result of the one-file cases."""
    if c["dt"] is None:
        return None
    dt = float.fromhex(c["dt"])
    result = c["result"] if len(c["files"]) == 1 else single.get((f, c["dt"]))
    if result is None and stored:  # times stored in the file do not depend on dt
        result = single.get((f, None))
    if result is None or len(result["frames"]) < 2 or "time" not in result["frames"][1]:
        return None
    step = float.fromhex(result["frames"][1]["time"]) - float.fromhex(result["frames"][0]["time"])
    if abs(step - dt) > 1e-4 * max(abs(dt), abs(step)):
        return f"inputs.dt_ps ({dt!r} ps) replaces the time step in {{path{k}}} ({step!r} ps)"
    return None


def trajectory_warnings(c, tmp, single):
    """The reader warnings (in order, "{pathK}" for the files) of a trajectory case MDBridge reads."""
    import scipy.io
    out = []
    for k, f in enumerate(c["files"]):
        path = tmp(f)
        p = "{path%d}" % k
        stored = False
        with open(path, "rb") as fh:
            head = fh.read(96)
        if head[:3] == b"CDF":
            if head[3] == 1:  # RD-N08
                out.append(p + ": NetCDF file uses 32-bit offsets (CDF-1); the Amber convention asks for "
                               "64-bit offsets (CDF-2)")
            try:
                with scipy.io.netcdf_file(path, "r", mmap=False) as nc:
                    frames = nc.variables["coordinates"].shape[0]
                    has_time = "time" in nc.variables
            except ValueError:  # a record count scipy cannot read (streaming.nc, RD-N05): the header's names
                with open(path, "rb") as fh:
                    has_time = b"\0\0\0\x04time" in fh.read(1 << 16)
                frames = None
            if len(c["files"]) == 1:
                frames = c["result"]["n_frames"]
            if not has_time and c["dt"] is None and frames > 1:  # RD-N10
                out.append(p + NO_TIMES)
            stored = has_time
        elif head[:4] in (b"\x54\0\0\0", b"\0\0\0\x54") and head[4:8] == b"CORD":
            nsavc = int.from_bytes(head[16:20], "little" if head[0] == 0x54 else "big", signed=True)
            if nsavc == 0:  # RD-D01
                out.append(p + ": the DCD header's frame interval (NSAVC) is 0; it is taken as 1 step")
        elif (c["format"] == "TRJ" or f.endswith((".mdcrd", ".trj"))) and len(c["files"]) == 1:  # RD-A04
            lines = open(path, encoding="utf-8", errors="replace").read().splitlines()
            while len(lines) > 1 and not lines[-1].strip():
                lines.pop()
            frames = c["result"]["frames"]
            per_frame = math.ceil(3 * c["atoms"] / 10) + (frames[0]["box"] is not None)
            first = 1 + len(frames) * per_frame
            if first < len(lines):
                out.append(f"{p}: the last frame is incomplete (lines {first + 1}-{len(lines)}) and is ignored")
        step = _time_step_warning(c, k, f, stored, single)
        if step:
            out.append(step)
    return out



# ------------------------------------------------------------------------- loader and export

def adjust_loader_summaries(summaries):
    """make_loader_cases.py's summaries (in place): RD-N02 and RD-A02 messages."""
    for case in summaries:
        error = case.get("error")
        if error is None:
            continue
        if error["message"] == ("Supplied n_atoms (7) != natom from ncdf (773). Note: n_atoms can be None and "
                                "then the ncdf value is used!"):
            error["message"] = ("The topology and NCDF trajectory files don't have the same number of atoms!\n"
                                f"Topology number of atoms 7\nTrajectory: {case['trajectories'][0]} Number of "
                                "atoms 773")
        elif error["message"] == "cannot reshape array of size 30 into shape (7,3)" and \
                case["trajectories"][0].endswith("ten.mdcrd"):
            # ten.mdcrd: 10 atoms, 3 lines of 10 per frame; 7 atoms need 10, 10, 1
            case["error"] = {"type": "OSError",
                             "message": f"{case['trajectories'][0]} line 4: expected 1 coordinate for 7 atoms, "
                                        "found 10 (does the topology match?)"}


_REPAIRED_DIR = None


def parm7_as_mdbridge_reads(path):
    """(file to load, fix(universe) or None): MDAnalysis' Universe of `path` as MDBridge reads it.

    chain_blank.parm7 (RD-P05): MDAnalysis drops the blank RESIDUE_CHAINID field of residue 3 and
    then ignores the section; MDBridge reads it as the chain "". The universe is loaded from a copy
    with the chain "Z", which the fix renames "" (chain IDs and segment IDs)."""
    global _REPAIRED_DIR
    if os.path.basename(path) != "chain_blank.parm7":
        return path, None
    import tempfile
    if _REPAIRED_DIR is None:
        import atexit
        import shutil
        _REPAIRED_DIR = tempfile.mkdtemp(prefix="mdbridge-repaired-")
        atexit.register(shutil.rmtree, _REPAIRED_DIR, True)
    with open(path, encoding="utf-8", newline="") as fh:
        text = fh.read()
    if "A   A       B   " not in text:
        raise RuntimeError("chain_blank.parm7: unexpected RESIDUE_CHAINID")
    copy = os.path.join(_REPAIRED_DIR, "chain_blank.parm7")
    with open(copy, "w", encoding="utf-8", newline="") as fh:
        fh.write(text.replace("A   A       B   ", "A   A   Z   B   "))

    def fix(u):
        u.atoms.chainIDs = np.array(["" if c == "Z" else c for c in u.atoms.chainIDs], dtype=object)
        u.segments.segids = np.array(["" if s == "Z" else s for s in u.segments.segids], dtype=object)

    return copy, fix
