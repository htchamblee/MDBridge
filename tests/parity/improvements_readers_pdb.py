"""MDBridge's improvements of the PDB / GRO readers and the loader (docs/IMPROVEMENTS.md 2.2, 2.3,
2.7; docs/DIFFERENCES.md "Readers and topologies"), applied to expected data on the test side.

The generators run MDAnalysis 2.11.0 / MDATK; where MDBridge deliberately does better, the
expected data is adjusted here, per item id:

* equivalent inputs: a file is rewritten into one that MDAnalysis reads the way MDBridge reads the
  original, and MDAnalysis is run on the rewrite (equivalent_file):
    RD-B07  bytes that are not UTF-8 are read as Latin-1 (re-encoded as UTF-8);
    RD-B04  lines that start with ATOM / HETATM but not in columns 1-6 ("ATOM  " / "HETATM") are
            ignored by both passes (dropped);
    RD-B03  CRYST1 records before the first MODEL count for the models without their own (moved
            into model 1);
    RD-B01 / RD-B09  CONECT records in the form MDBridge reads them (hybrid-36 serials, a blank
            partner field is no partner, a malformed record is skipped alone): rewritten in plain
            decimal where they fit, dropped where not (the dumps' bonds come from conect_bonds());
    RD-B10  atoms with a blank element column in a file whose column is partly filled get the
            element guessed from their name (with MDBridge's residue rule), as in a file without
            the column;
    RD-G01 / RD-G05 / RD-G07  GRO: coordinates re-aligned where MDAnalysis' width rule fails, box
            values as MDBridge reads them, one file per frame;
* patched results: conect_bonds() (RD-B01, RD-B09), the PDB reader's HEADER / COMPND / REMARK
  texts sliced by characters of the decoded line (RD-B07; MDAnalysis slices bytes), the cells of
  frames read out of order (RD-B03: a frame's cell no longer depends on the frame read before),
  every frame of a multi-frame GRO file (RD-G07);
* errors MDBridge words under the reader-message rule (ERRORS: RD-B02, RD-B05, RD-B06, RD-B08,
  RD-G02, RD-G03, RD-G04, RD-G05) and the loader's (LD-04, LD-06), chains of mixed formats read
  per file (LD-01);
* the command line (cli_cases.json): the reader warnings `info` / `export-pdb` / `analyze` now print
  (LD-03, LD-05) and the clearer errors (LD-06).

Pure Python (MDAnalysis only through the callbacks and Universe objects the generators pass in).
"""

import codecs
import os
import random
import re

# --- text ------------------------------------------------------------------------------------------


def _latin1_fallback(exc):
    """RD-B07: a byte that is not part of valid UTF-8 is that Latin-1 character."""
    return chr(exc.object[exc.start]), exc.start + 1


codecs.register_error("mdbridge_latin1", _latin1_fallback)


def decode(raw: bytes) -> str:
    return raw.decode("utf-8", errors="mdbridge_latin1")


LINE_RE = re.compile(r"[^\r\n]*(?:\r\n|\r|\n)|[^\r\n]+$")


def lines_of(text: str):
    """The lines of a text with their terminators (\\n, \\r\\n, \\r end lines)."""
    return LINE_RE.findall(text)


def body(line: str) -> str:
    return line.rstrip("\r\n")


def is_atom_record(line: str) -> bool:
    return line[:6] in ("ATOM  ", "HETATM")


# --- PDB: CONECT as MDBridge reads it (RD-B01, RD-B09) ---------------------------------------------

def _hy36(field: str):
    """Columns 7-11 as a serial: an int, else hybrid-36 of 5 characters; None when not a serial."""
    try:
        return int(field)
    except ValueError:
        pass
    if len(field) != 5 or not ("A" <= field[0] <= "Z" or "a" <= field[0] <= "z"):
        return None
    upper = field[0] <= "Z"
    digits = "0123456789" + ("ABCDEFGHIJKLMNOPQRSTUVWXYZ" if upper else "abcdefghijklmnopqrstuvwxyz")
    value = 0
    for c in field:
        if c not in digits:
            return None
        value = value * 36 + digits.index(c)
    return value - 10 * 36 ** 4 + 100000 if upper else value + 16 * 36 ** 4 + 100000


def conect_records(lines):
    """[(base, [partners])] of the CONECT records MDBridge reads; malformed ones are left out."""
    out = []
    for raw in lines:
        if raw[:6] != "CONECT":
            continue
        line = raw.strip()
        base = _hy36(line[6:11])
        length = max(len(line) - 11, 0)
        count = length // 5
        if base is None or count == 0 or length % count != 0:
            out.append((base, None if base is None or (count and length % count) else []))
            continue
        partners = []
        for k in range(count):
            field = line[11 + 5 * k:16 + 5 * k]
            if not field.strip():
                continue
            value = _hy36(field)
            if value is None:
                partners = None
                break
            partners.append(value)
        out.append((base, partners))
    return [(b, p) for b, p in out if b is not None and p is not None]


def conect_bonds(lines, serials):
    """The bonds MDBridge reads: sorted [i, j] (i < j), no self-bonds; None without any pair."""
    atom_of = {s: i for i, s in enumerate(serials)}
    bonds, any_pair = set(), False
    for base, partners in conect_records(lines):
        for p in partners:
            if base in atom_of and p in atom_of:
                any_pair = True
                i, j = atom_of[base], atom_of[p]
                if i != j:
                    bonds.add((min(i, j), max(i, j)))
    return sorted(map(list, bonds)) if any_pair else None


# --- PDB: equivalent inputs ------------------------------------------------------------------------

def topology_atom_lines(lines):
    """Indices of the lines the topology pass reads as atoms (up to the first END line)."""
    out = []
    for k, line in enumerate(lines):
        s = line.strip()
        if not s:
            continue
        if s.startswith("END"):
            break
        if is_atom_record(line):
            out.append(k)
    return out


def _canonical_conect(lines):
    """RD-B01 / RD-B09: CONECT records as plain decimal records MDAnalysis parses as MDBridge does."""
    out = []
    for line in lines:
        if line[:6] != "CONECT":
            out.append(line)
            continue
        single = conect_records([line])
        if not single:
            continue  # malformed: skipped
        base, partners = single[0]
        if not all(0 <= v <= 99999 for v in [base, *partners]):
            continue  # not expressible in decimal; the dumps take conect_bonds()
        ending = line[len(body(line)):] or "\n"
        out.append("CONECT" + "".join(f"{v:5d}" for v in [base, *partners]) + ending)
    return out


def _move_preamble_cryst(lines):
    """RD-B03: a CRYST1 before the first MODEL (when the CRYST1 and MODEL counts differ) goes into
    model 1 when that model has none of its own."""
    models = [k for k, line in enumerate(lines) if line.startswith("MODEL")]
    crysts = [k for k, line in enumerate(lines) if line.startswith("CRYST")]
    if not models or len(models) == len(crysts):
        return lines
    preamble = [k for k in crysts if k < models[0] and lines[k][:6] == "CRYST1"]
    end = models[1] if len(models) > 1 else len(lines)
    if not preamble or any(models[0] < k < end and lines[k][:6] == "CRYST1" for k in crysts):
        return lines
    last = preamble[-1]
    moved = lines[last]
    out = [line for k, line in enumerate(lines) if k != last]
    at = models[0]  # its index after the removal
    out.insert(at, moved)
    return out


RECORD_COLUMNS = {"HEADER": (10, 66), "TITLE": (8, 80), "COMPND": (7, 80), "REMARK": (6, None)}


def _slices_agree(line: str) -> bool:
    """Whether MDAnalysis' byte slice of a HEADER / TITLE / COMPND / REMARK record decodes to
    MDBridge's character slice."""
    for record, (first, last) in RECORD_COLUMNS.items():
        if line.startswith(record):
            raw = line.encode("utf-8")
            try:
                return raw[first:last].decode("utf-8").strip() == line[first:last].strip()
            except UnicodeDecodeError:
                return False
    return True


def equivalent_pdb_text(raw: bytes) -> list:
    """The original's lines after the text-level rewrites (RD-B07, RD-B04, RD-B03, CONECT)."""
    lines = lines_of(decode(raw))
    lines = [line for line in lines
             if is_atom_record(line) or not line.strip().startswith(("ATOM", "HETATM"))]
    lines = _move_preamble_cryst(lines)
    # MDAnalysis decodes byte slices of these records (a cut character fails); MDBridge slices the
    # decoded line: patch_reader() takes their texts from the original
    lines = [line if _slices_agree(line) else re.sub(r"[^\x00-\x7f]", "?", line) for line in lines]
    return _canonical_conect(lines)


def _fill_elements(lines, guess_types):
    """RD-B10: blank element columns of a partly filled column get the guessed element."""
    atoms = topology_atom_lines(lines)
    columns = [body(lines[k])[76:78].strip() for k in atoms]
    if not any(columns) or all(columns):
        return lines
    blanked = list(lines)
    for k in atoms:
        text = body(lines[k])
        blanked[k] = (text[:76].ljust(76) + "  " + text[78:] if len(text) > 76 else text) + \
            lines[k][len(text):]
    types = guess_types(blanked)
    out = list(lines)
    for k, column, guessed in zip(atoms, columns, types):
        if column or len(guessed) > 2:
            continue
        text = body(lines[k])
        out[k] = text[:76].ljust(76) + guessed.rjust(2) + text[78:] + lines[k][len(text):]
    return out


# --- GRO: equivalent inputs ------------------------------------------------------------------------

def gro_box_values(line):
    """MDBridge's box line (RD-G05): split; else 10-character fields; else run-together numbers
    with five decimals (signs kept); None when none gives 3 or 9 values."""
    text = line.rstrip()

    def numbers(tokens):
        try:
            values = [float(t) for t in tokens]
        except ValueError:
            return None
        return values if len(values) in (3, 9) else None

    for tokens in (text.split(),
                   [text[i:i + 10] for i in range(0, len(text), 10)] if text and len(text) % 10 == 0 else ["x"],
                   re.fullmatch(r"(?:\s*-?\d+\.\d{5})+\s*", text) and re.findall(r"-?\d+\.\d{5}", text) or ["x"]):
        values = numbers(tokens)
        if values:
            return tokens
    return None


def gro_frames(raw: bytes):
    """The frames of a GRO file (RD-G07) as lists of lines with terminators."""
    lines = lines_of(decode(raw))
    n = int(lines[1])
    frames, at = [], 0
    while at < len(lines) and "".join(lines[at:]).strip():
        block = lines[at:at + n + 3]
        if len(block) < n + 3:
            break
        try:
            if int(block[1]) != n:
                break
        except ValueError:
            break
        frames.append(block)
        at += n + 3
    return frames


def _realign_gro(block):
    """RD-G01 / RD-G05: where MDAnalysis' width rule (the first '.' from column 26) misses
    MDBridge's width (the distance between the decimal points of y and z), the atom fields are
    rewritten with the decimal point 4 characters into each field; the box space-separated."""
    first = body(block[2])
    dots = [20 + m.start() for m in re.finditer(r"\.", first[20:])][:3]
    width = dots[2] - dots[1]
    found = first.find(".", 25)
    out = list(block[:2])
    for line in block[2:-1]:
        text = body(line)
        if found - 24 == width:
            out.append(line)
            continue
        fields = [text[20 + k * width:20 + (k + 1) * width].strip() for k in range(3)]
        fields = [f.split(".")[0].rjust(4) + "." + f.split(".")[1] for f in fields]
        out.append(text[:20] + "".join(fields) + "\n")
    box = gro_box_values(body(block[-1]))
    out.append(" ".join(box) + "\n" if box else block[-1])
    return out


# --- files -----------------------------------------------------------------------------------------

def guess_types(lines):
    """The types MDAnalysis (with MDBridge's tables and residue rule, independent_tables.py)
    guesses for the atoms of a PDB text (a file without element column)."""
    import tempfile
    import MDAnalysis as mda
    import independent_tables
    independent_tables.apply()
    with tempfile.TemporaryDirectory(prefix="guess-") as tmp:
        path = os.path.join(tmp, "guess.pdb")
        with open(path, "w", encoding="utf-8", newline="") as fh:
            fh.write("".join(lines))
        u = mda.Universe(path)
        independent_tables.residue_rule_universe(u, path)
        return [str(t) for t in u.atoms.types]


def equivalent_file(path: str, folder: str, guess=guess_types) -> str:
    """A file in `folder` with the same name that MDAnalysis reads as MDBridge reads `path`.
    guess(lines) -> the types guessed for the atoms of a PDB text without element column."""
    with open(path, "rb") as fh:
        raw = fh.read()
    target = os.path.join(folder, os.path.basename(path))
    if path.lower().endswith(".gro"):
        frames = gro_frames(raw) if raw else []
        text = "".join(_realign_gro(frames[0])) if frames else decode(raw)
    elif raw:
        lines = equivalent_pdb_text(raw)
        lines = _fill_elements(lines, guess)
        text = "".join(lines)
    else:
        text = ""
    with open(target, "w", encoding="utf-8", newline="") as fh:
        fh.write(text)
    return target


def gro_frame_files(path: str, folder: str) -> list:
    """One equivalent single-frame GRO file per frame (RD-G07)."""
    with open(path, "rb") as fh:
        frames = gro_frames(fh.read())
    out = []
    for k, block in enumerate(frames):
        sub = os.path.join(folder, f"frame{k}")
        os.makedirs(sub, exist_ok=True)
        target = os.path.join(sub, os.path.basename(path))
        with open(target, "w", encoding="utf-8", newline="") as fh:
            fh.write("".join(_realign_gro(block)))
        out.append(target)
    return out


def patch_dump(data: dict, path: str) -> dict:
    """RD-B01 / RD-B09: the bonds of a PDB dump as MDBridge reads the CONECT records."""
    if not path.lower().endswith(".pdb"):
        return data
    with open(path, "rb") as fh:
        lines = lines_of(decode(fh.read()))
    serials = data["atoms"]["ids"]
    made_up = any(_hy36(body(lines[k]).strip()[6:11]) is None for k in topology_atom_lines(lines))
    data["bonds"] = None if made_up else conect_bonds(lines, serials)
    return data


def frame_order(n):
    """The frames reader_info reads: in order, then a shuffled order (seed 3)."""
    order = list(range(n))
    shuffled = order[::-1] + order
    random.Random(3).shuffle(shuffled)
    return order + shuffled


def patch_reader(info: dict, path: str, frames_of=None) -> dict:
    """The PDB / GRO reader's results as MDBridge gives them: records sliced by characters
    (RD-B07), cells independent of the read order (RD-B03), every GRO frame (RD-G07:
    frames_of(dt) -> [frame dicts of frame 0 of each frame file])."""
    if path.lower().endswith(".gro"):
        if frames_of is not None:
            per_frame = frames_of(info["dt_given"])
            if len(per_frame) > 1:
                dt = info["dt_given"] or 1.0
                info["n_frames"] = len(per_frame)
                info["frames"] = [dict(per_frame[k], frame=k, time=float(k * dt).hex())
                                  for k in frame_order(len(per_frame))]
        return info
    with open(path, "rb") as fh:
        raw = fh.read()
    header, compounds, remarks = "", [], []
    for line in raw.split(b"\n"):
        text = decode(line + b"\n")
        if line.startswith(b"HEADER"):
            header = text[10:66].strip()
        elif line.startswith(b"COMPND"):
            compounds.append(text[7:80].strip())
        elif line.startswith(b"REMARK"):
            remarks.append(text[6:].strip())
    info["header"], info["compound"], info["remarks"] = header, compounds, remarks
    frames = info["frames"]
    in_order = {f["frame"]: f for f in frames[:info["n_frames"]] if "error" not in f}
    for f in frames[info["n_frames"]:]:
        if "error" not in f and f["frame"] in in_order:
            f["dimensions"] = in_order[f["frame"]]["dimensions"]
    return info


# --- errors ----------------------------------------------------------------------------------------

def _error(message, kind="ValueError"):
    return {"type": kind, "message": message}


# file -> the error of reading it (topology and its coordinates), "{p}" for its path as given
ERRORS = {
    "hy36_bad.pdb": "{p} line 1: invalid hybrid-36 atom serial 'Aa000'",                  # RD-B02
    "empty.pdb": "{p}: file is empty",                                                    # RD-B06
    "remarks_only.pdb": "{p}: no ATOM or HETATM records",                                 # RD-B05
    "endmdl_first.pdb": "{p} line 5: ATOM or HETATM record after END (line 4)",           # RD-B05
    "models_mismatch.pdb": "{p}: model 2 has 6 atoms, model 1 has 7 (line 11)",           # RD-B05
    "models_more.pdb": "{p}: model 2 has 6 atoms, model 1 has 5 (line 9)",                # RD-B05
    "bad_coordinate.pdb": "{p} line 3: x coordinate '   1.0x0' is not a number",          # RD-B08
    "no_box.gro": "{p}: missing box line after atom 7 (line 9)",                          # RD-G02
    "too_few.gro": "{p}: the header (line 2) says 7 atoms but the file has 4 atom lines",  # RD-G03
    "bad_count.gro": "{p} line 2: atom count ' seven' is not an integer",                 # RD-G04
    "bad_box.gro": "{p} line 10: GRO box line '   3.0   4.0' has neither 3 nor 9 values",  # RD-G05
}


def expected_error(name: str, path_text: str):
    """The error MDBridge gives for the file (path_text: its path as the messages show it)."""
    message = ERRORS.get(name)
    return None if message is None else _error(message.replace("{p}", path_text))


def loader_error(topology: str, trajectories, dt):
    """LD-04 / LD-06: errors of load_universe MDBridge words itself (None: as MDATK)."""
    if dt is not None and dt == dt and dt >= 0 and dt != float("inf"):
        pass
    elif dt is not None:
        return _error("inputs.dt_ps must be a positive number of picoseconds (0: the files' time step), got "
                      + repr(float(dt)))
    if not trajectories and topology.lower().endswith((".parm7", ".prmtop", ".top")):
        return _error(f"{topology} has no coordinates; give a trajectory with --traj")
    return None


MDATK_FORMATS = {".nc": "NCDF", ".ncdf": "NCDF", ".netcdf": "NCDF", ".mdcrd": "TRJ", ".crd": "TRJ", ".trj": "TRJ",
                 ".rst7": "INPCRD", ".inpcrd": "INPCRD", ".ncrst": "NCRST"}


def chain_with_formats(trajectories):
    """LD-01: a chain whose files have different extensions, each with MDATK's format for its
    extension (ChainReader's (file, format) pairs); None for chains MDATK already reads right."""
    exts = {os.path.splitext(t)[1].lower() for t in trajectories}
    if len(trajectories) < 2 or len(exts) < 2 or not all(os.path.exists(t) for t in trajectories):
        return None
    out = []
    for t in trajectories:
        fmt = MDATK_FORMATS.get(os.path.splitext(t)[1].lower())
        out.append((t, fmt) if fmt else t)
    return out


def time_step_warnings(trajectories, dt):
    """LD-05: the warning MDBridge gives for a file that keeps its own times although inputs.dt_ps is
    given (NetCDF with a time variable; every other reader then takes frame * dt)."""
    out = []
    if not dt or dt != dt:
        return out
    import MDAnalysis as mda  # the generators' MDAnalysis
    for path in trajectories:
        if os.path.splitext(path)[1].lower() not in (".nc", ".ncdf", ".netcdf"):
            continue
        reader = mda.coordinates.TRJ.NCDFReader(path)
        try:
            if reader.n_frames < 2 or "time" not in reader.trjfile.variables:
                continue
            times = reader.trjfile.variables["time"]
            step = float(times[1]) - float(times[0])
        finally:
            reader.close()
        if abs(step - dt) > 1e-4 * max(abs(dt), abs(step)):
            out.append(f"inputs.dt_ps ({float(dt)!r} ps) differs from the time step in {path} ({step!r} ps); "
                       "the file's times are used")
    return out


# --- the command line (cli_cases.json) ---------------------------------------------------------------

DT_WARNING = ("WARNING: inputs.dt_ps ({dt} ps) differs from the time step in {{work}}/demo_dry.nc (10.0 ps); "
              "the file's times are used\n")

ELEMENTS_WARNINGS = (
    "WARNING: {work}/elements.pdb: 2 atoms have no element symbol in columns 77-78; guessed from the atom name "
    "(first: line 6)\n"
    "WARNING: {work}/elements.pdb: 3 atoms have an unknown element symbol in columns 77-78; the element is left "
    "blank (first: 'XX' on line 4)\n"
    "WARNING: {work}/elements.pdb: 3 atoms have no known element mass and get mass 0 (first: 'X1' on line 4)\n")


CHAINS_WARNINGS = (
    "WARNING: {work}/chains.pdb line 6: formal charge '-' in columns 79-80 is not a charge (such as '2-'); no "
    "formal charges are read\n"
    "WARNING: {work}/chains.pdb: 2 atoms have an unknown element symbol in columns 77-78; the element is left "
    "blank (first: 'l1' on line 6)\n"
    "WARNING: {work}/chains.pdb: 2 atoms have no known element mass and get mass 0 (first: 'CL' on line 6)\n")


def adjust_cli(name: str, argv, res: dict) -> dict:
    """The expected result of a command-line case with MDBridge's reader warnings and errors."""
    stderr = res.get("stderr")
    dt = None
    for k, a in enumerate(argv):
        if a == "-s" and argv[k + 1].startswith("inputs.dt_ps="):
            dt = float(argv[k + 1].split("=", 1)[1])
    uses_nc = "demo_dry.nc" in argv or any(a in ("-p", "--project") for a in argv)
    if dt is None and any(a in ("-p", "--project") for a in argv):
        dt = 5.0  # the cases' project file (PROJECT in cli_parity.py)
    if name == "info_parm7_alone":  # LD-06
        res["stderr"] = "error: {work}/demo_dry.parm7 has no coordinates; give a trajectory with --traj\n"
        res.pop("mdatk_crash", None)
        return res
    # the commands that load a system; analyze and export-pdb log through their reporter (-q: silent)
    loads = argv[0] == "info" or argv[0] in ("analyze", "export-pdb") and not ({"-q", "--quiet"} & set(argv))
    if stderr is None or not loads:
        return res
    lines = stderr.splitlines(keepends=True)
    warnings = ""
    if dt and uses_nc and dt != 10.0:  # LD-05
        warnings += DT_WARNING.format(dt=repr(dt))
    if name == "info_pdb":  # LD-03, RD-B10, RD-B11
        warnings += ELEMENTS_WARNINGS
    if name == "export_pdb_trajectory":  # LD-03, RD-B11 (io/pdb/export/files/chains.pdb)
        warnings += CHAINS_WARNINGS
    if not warnings:
        return res
    # the readers' warnings come when the system is loaded: after "inputs ready", else first
    at = next((k + 1 for k, line in enumerate(lines) if line.endswith("inputs ready\n")), 0)
    res["stderr"] = "".join(lines[:at]) + warnings + "".join(lines[at:])
    return res
