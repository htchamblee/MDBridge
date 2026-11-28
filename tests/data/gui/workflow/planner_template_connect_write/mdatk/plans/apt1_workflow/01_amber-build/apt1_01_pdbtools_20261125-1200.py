#!/usr/bin/env python3
"""PDB hygiene for D- and L-RNA structures - copied by MDATK <VERSION>.

Usage: python3 THIS_SCRIPT clean IN.pdb -o OUT.pdb --scheme amber [--l-ranges C:45-73]
       python3 THIS_SCRIPT check IN.pdb
"""

import copy
import math
import os
import re
import string
import sys
from collections import Counter
from dataclasses import dataclass, field
from typing import Dict, Iterator, List, NamedTuple, Optional, Sequence, Tuple, Union

Vec = Tuple[float, float, float]


class PDBError(ValueError):
    """Unreadable input or an impossible request (e.g. a chain-restore mismatch)."""


CHAIN_IDS = string.ascii_uppercase + string.ascii_lowercase + string.digits
SCHEMES = ("amber", "ccd", "rosetta", "pdb")
RECORD_POLICIES = ("keep", "rosetta", "pdb")

# Signs of the chiral volumes of an L-ribose: (centre, (atom 1, atom 2, atom 3), sign), from the
# _chem_comp_chir tables of the 0A/0C/0G/0U restraint dictionaries; "N" is N9 (purines) or N1
# (pyrimidines).  V = a.(b x c) with a, b, c the vectors from the centre to atoms 1-3; D has the
# opposite signs.
MIN_CHIRAL_VOLUME = 0.5  # Å^3; real sugar centres give about 1.8-2.6, schematic models can be near 0
L_CHIRAL_SIGNS = (("C4'", ("O4'", "C3'", "C5'"), 1.0),
                  ("C3'", ("O3'", "C4'", "C2'"), -1.0),
                  ("C1'", ("O4'", "N", "C2'"), 1.0))

_AMINO = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G",
          "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S",
          "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"}
# Amber protonation / disulfide variants: one-letter codes only, they are not standard PDB names
_AMINO_VARIANTS = {"HID": "H", "HIE": "H", "HIP": "H", "CYX": "C", "CYM": "C", "ASH": "D", "GLH": "E",
                   "LYN": "K"}
_STANDARD_NAMES = frozenset(_AMINO) | {"A", "C", "G", "U", "I", "N", "DA", "DC", "DG", "DT", "DI", "DN",
                                       "UNK"}
_WATER = frozenset({"HOH", "WAT", "H2O", "DOD", "SOL", "TIP", "TIP3", "TIP4", "TIP5", "T3P", "T4P", "T4E",
                    "SPC", "OPC"})
_NT_RE = re.compile(r"^([RL0D]?)([ACGTU])([53N]?)$")
_NT_ALIASES = {"ADE": ("A", "rna"), "CYT": ("C", "rna"), "GUA": ("G", "rna"), "URA": ("U", "rna"),
               "URI": ("U", "rna"), "THY": ("T", "dna")}
_TWO_LETTER_ELEMENTS = frozenset({"MG", "NA", "CL", "ZN", "MN", "CA", "LI", "RB", "CS", "SR", "BA", "CD",
                                  "FE", "CO", "NI", "CU", "BR", "AG", "HG", "PB", "TL", "AU", "PT"})
_SEQRES_ORDER = ("MODRES", "HET", "HETNAM", "HETSYN", "FORMUL", "HELIX", "SHEET", "SSBOND", "LINK",
                 "CISPEP", "SITE", "CRYST1", "ORIGX", "SCALE", "MTRIX")
_COORD_RE = re.compile(r"-?\d+\.\d{3}")
_FLOAT_RE = re.compile(r"-?\d+\.\d+")


# --- data model ---------------------------------------------------------------------------------------

@dataclass
class Atom:
    """One ATOM/HETATM record; ``tail`` is columns 28-76 verbatim (coordinates, occupancy, B, segid)."""
    record: str
    name: str
    alt: str
    tail: str
    element: str = ""
    charge: str = ""
    serial: Optional[int] = None  # input serial, used to remap kept CONECT records
    anisou: str = ""  # columns 28-80 of a following ANISOU record

    @property
    def xyz(self) -> Vec:
        t = self.tail
        return (float(t[3:11]), float(t[11:19]), float(t[19:27]))


@dataclass
class Residue:
    name: str
    chain: str
    resseq: int
    icode: str = ""
    atoms: List[Atom] = field(default_factory=list)
    chirality: Optional[str] = None  # "D" / "L" once assign_chirality ran (override, detection or name)

    @property
    def label(self) -> str:
        return f"{self.chain or '_'}:{self.resseq}{self.icode} {self.name}"

    def coords(self) -> Dict[str, Vec]:
        """Atom name -> coordinates (first alternate location; old-style '*' read as a prime)."""
        out: Dict[str, Vec] = {}
        for a in self.atoms:
            out.setdefault(a.name.replace("*", "'"), a.xyz)
        return out


@dataclass
class Structure:
    """A parsed PDB file: header lines, TER-delimited segments of residues, trailing lines."""
    header: List[str] = field(default_factory=list)
    segments: List[List[Residue]] = field(default_factory=list)
    trailer: List[str] = field(default_factory=list)
    conect: List[str] = field(default_factory=list)
    model: str = ""  # MODEL line of the kept model ("" when the input had none)
    notes: List[str] = field(default_factory=list)  # parse notes

    def residues(self) -> Iterator[Residue]:
        for seg in self.segments:
            yield from seg

    def chain_ids(self) -> List[str]:
        """Chain ids in order of first appearance."""
        return list(dict.fromkeys(r.chain for r in self.residues()))


# --- residue names ------------------------------------------------------------------------------------

class ResidueName(NamedTuple):
    base: str  # A C G U (T for DNA)
    sugar: str  # "rna" / "dna"
    chirality: str  # "D" / "L" as implied by the name
    terminal: str  # "5", "3", "N" (single residue) or ""


def parse_residue_name(name: str) -> Optional[ResidueName]:
    """Recognise a nucleotide name in any supported scheme or common alias; None otherwise.

    RNA: A, A5, A3, AN, RA, RA5, RA3, RAN, ADE/CYT/GUA/URA/URI; L-RNA: LA, LA5, LA3, 0A, 0A5, 0A3 (and the
    N forms); DNA: DA, DA5, DA3, DAN, T, THY ...
    """
    n = name.strip().upper()
    if n in _NT_ALIASES:
        base, sugar = _NT_ALIASES[n]
        return ResidueName(base, sugar, "D", "")
    m = _NT_RE.match(n)
    if not m:
        return None
    pre, base, term = m.groups()
    if pre == "D" or (base == "T" and pre == ""):
        return ResidueName(base, "dna", "D", term)
    if base == "T":
        return None
    return ResidueName(base, "rna", "L" if pre in ("L", "0") else "D", term)


def scheme_name(base: str, chirality: str, terminal: str, scheme: str) -> str:
    """Residue name of an RNA nucleotide in ``scheme`` (terminal: "5", "3", "N" or "")."""
    if scheme == "amber":
        if chirality == "L":  # nucleic12_LRNA.lib defines LX, LX5 and LX3 only
            return "L" + base + {"5": "5", "3": "3", "N": "5"}.get(terminal, "")
        return base + terminal
    if scheme == "ccd":
        return "0" + base + {"5": "5", "3": "3"}.get(terminal, "") if chirality == "L" else base
    if scheme in ("rosetta", "pdb"):
        return "0" + base if chirality == "L" else base
    raise PDBError(f"unknown residue name scheme {scheme!r} (choose from {', '.join(SCHEMES)})")


def _nucleotide(res: Residue) -> Optional[ResidueName]:
    """Name-based recognition, guarded by backbone atoms so that e.g. a lanthanum ion 'LA' is not RNA."""
    info = parse_residue_name(res.name)
    if info is None or not any("'" in a.name or "*" in a.name or a.name == "P" for a in res.atoms):
        return None
    return info


def _is_polymer(res: Residue) -> bool:
    return _nucleotide(res) is not None or res.name.upper() in _AMINO or res.name.upper() in _AMINO_VARIANTS


def _is_solvent(res: Residue) -> bool:
    """Water or a single-atom residue (ion)."""
    return not _is_polymer(res) and (res.name.upper() in _WATER or len(res.atoms) == 1)


def _code(res: Residue) -> str:
    """One-letter code: nucleotide base, amino acid, else 'X'."""
    nt = _nucleotide(res)
    if nt:
        return nt.base
    return _AMINO.get(res.name.upper()) or _AMINO_VARIANTS.get(res.name.upper()) or "X"


def _chirality(res: Residue) -> Optional[str]:
    nt = _nucleotide(res)
    return res.chirality or (nt.chirality if nt else None)


# --- parsing ------------------------------------------------------------------------------------------

def _read_text(source) -> str:
    if isinstance(source, os.PathLike) or (isinstance(source, str) and "\n" not in source
                                           and os.path.isfile(source)):
        with open(source) as fh:
            return fh.read()
    if not isinstance(source, str):
        raise PDBError(f"expected PDB text or a path, got {type(source).__name__}")
    return source


def _fused_field(line: str) -> Optional[Tuple[str, str, str, int]]:
    """ChimeraX 'resname+chain:' field -> (alt, resname, chain, colon index), or None for a normal line."""
    c = line.find(":", 18, 30)
    if c < 0 or line[c - 1].isspace():
        return None
    words = line[16:c - 1].split()  # [resname] or [altloc, resname]
    if len(words) == 1:
        return " ", words[0], line[c - 1], c
    if len(words) == 2 and words[0] == line[16]:
        return words[0], words[1], line[c - 1], c
    return None


def _tokens_tail(rest: str, lineno: int) -> Tuple[int, str, str, str]:
    """Fallback for shifted columns: resseq(+icode), x, y, z, occupancy, B, element from the tokens."""
    m = re.match(r"\s*(-?\d+)([A-Za-z]?)\s", rest)
    coords = list(_COORD_RE.finditer(rest, m.end())) if m else []
    if len(coords) < 3:
        raise PDBError(f"line {lineno}: cannot read the damaged ATOM record")
    x, y, z = (float(c.group()) for c in coords[:3])
    after = rest[coords[2].end():]
    nums = [float(v) for v in _FLOAT_RE.findall(after)[:2]]
    occ, bfac = (nums + [1.0, 0.0][len(nums):])[:2]
    words = after.split()
    element = words[-1] if words and words[-1].isalpha() and len(words[-1]) <= 2 else ""
    tail = f"   {x:8.3f}{y:8.3f}{z:8.3f}{occ:6.2f}{bfac:6.2f}".ljust(49)
    return int(m.group(1)), m.group(2), tail, element


def _coords_ok(tail: str) -> bool:
    try:
        float(tail[3:11]), float(tail[11:19]), float(tail[19:27])
    except ValueError:
        return False
    return True


def _parse_atom(line: str, lineno: int) -> Tuple[Atom, str, str, int, str, bool]:
    line = line.ljust(80)
    fused = _fused_field(line)
    element = charge = ""
    if fused:
        alt, resname, chain, colon = fused
        rest = line[colon + 1:].ljust(58)  # aligned like columns 23-80 when the shift is uniform
        try:
            resseq, icode, tail = int(rest[0:4]), rest[4].strip(), rest[5:54]
            element, charge = rest[54:56].strip(), rest[56:58].strip()
            if not _coords_ok(tail):
                raise ValueError
        except ValueError:
            resseq, icode, tail, element = _tokens_tail(rest, lineno)
            charge = ""
    else:
        alt, resname, chain = line[16], line[17:21].strip(), line[21]
        try:
            resseq = int(line[22:26])
        except ValueError:
            raise PDBError(f"line {lineno}: unreadable residue number {line[22:26]!r}") from None
        icode, tail = line[26].strip(), line[27:76]
        element, charge = line[76:78].strip(), line[78:80].strip()
        if not _coords_ok(tail):
            raise PDBError(f"line {lineno}: unreadable coordinates")
    try:
        serial: Optional[int] = int(line[6:11])
    except ValueError:
        serial = None
    atom = Atom(record=line[:6].strip(), name=line[12:16].strip(), alt=alt.strip(), tail=tail,
                element=element, charge=charge, serial=serial)
    return atom, resname, chain.strip(), resseq, icode, bool(fused)


def read_pdb(source) -> Structure:
    """Parse PDB text or a path (ChimeraX damage repaired, first model only); notes in ``.notes``."""
    if isinstance(source, Structure):
        return copy.deepcopy(source)
    text = _read_text(source)
    s = Structure()
    seg: List[Residue] = []
    res: Optional[Residue] = None
    started = ended = skipping = False
    n_models = repaired = dropped = after_end = 0

    def close_segment():
        nonlocal seg, res
        if seg:
            s.segments.append(seg)
        seg, res = [], None

    for lineno, line in enumerate(text.splitlines(), 1):
        line = line.rstrip()
        rec = line[:6].strip().upper()
        if not line.strip():
            continue
        if ended:
            after_end += rec in ("ATOM", "HETATM")
            continue
        if rec == "MODEL":
            n_models += 1
            if n_models == 1:
                s.model = line
            else:
                skipping = True
            continue
        if rec == "ENDMDL":
            close_segment()
            skipping = True
            continue
        if rec in ("ATOM", "HETATM", "TER", "ANISOU", "SIGATM", "SIGUIJ") and skipping:
            continue
        if rec in ("ATOM", "HETATM"):
            started = True
            atom, resname, chain, resseq, icode, fixed = _parse_atom(line, lineno)
            repaired += fixed
            if seg and seg[-1].chain != chain:
                close_segment()
            if res is None or (res.name, res.chain, res.resseq, res.icode) != (resname, chain, resseq, icode):
                res = Residue(resname, chain, resseq, icode)
                seg.append(res)
            res.atoms.append(atom)
        elif rec == "TER":
            close_segment()
        elif rec == "ANISOU" and res is not None:
            res.atoms[-1].anisou = line.ljust(80)[27:80]
        elif rec in ("SIGATM", "SIGUIJ", "ANISOU"):
            dropped += 1
        elif rec == "CONECT":
            s.conect.append(line)
        elif rec == "MASTER":
            s.notes.append("Dropped the MASTER record (its counts no longer apply)")
        elif rec == "END":
            close_segment()
            ended = True
        elif started:
            s.trailer.append(line)
        else:
            s.header.append(line)
    close_segment()
    if not s.segments:
        raise PDBError("no ATOM/HETATM records found (not a PDB file, or the path does not exist)")
    if repaired:
        s.notes.append(f"Repaired {repaired} ChimeraX-damaged ATOM/HETATM lines (residue name fused with "
                       f"the chain id and a colon, e.g. 'GC:' -> residue G, chain C)")
    if n_models > 1:
        s.notes.append(f"Kept model 1 of {n_models}")
    if after_end:
        s.notes.append(f"Ignored {after_end} ATOM/HETATM lines after END (only the first model is kept)")
    if dropped:
        s.notes.append(f"Dropped {dropped} SIGATM/SIGUIJ/orphan ANISOU records")
    return s


# --- writing ------------------------------------------------------------------------------------------

def _element(atom: Atom, resname: str) -> str:
    if atom.element:
        return atom.element
    name = atom.name.strip()
    bare = re.sub(r"[\d+\-]", "", name).upper()
    if bare in _TWO_LETTER_ELEMENTS and bare == re.sub(r"[\d+\-]", "", resname).upper():
        return bare  # ion residue named after its element (MG, Na+, Cl- ...)
    letters = re.sub(r"^\d+", "", name)
    return letters[:1].upper() if letters else ""


def _atom_name_field(name: str, element: str) -> str:
    """PDB alignment: 4-character names and two-letter elements start in column 13, others in 14."""
    if len(name) >= 4 or len(element) == 2:
        return name[:4].ljust(4)
    return " " + name.ljust(3)


def _residue_field(resname: str) -> str:
    """Columns 18-21: names of up to 3 characters right-justified in 18-20, 4-character names in 18-21."""
    return f"{resname:>3} " if len(resname) <= 3 else resname[:4]


def _id_columns(res: Residue) -> str:
    """Columns 18-27 (residue name, chain, residue number, insertion code); numbers above 9999 wrap."""
    resseq = res.resseq % 10000 if res.resseq > 9999 else res.resseq
    return f"{_residue_field(res.name)}{(res.chain or ' ')[:1]}{resseq:>4}{(res.icode or ' ')[:1]}"


def _atom_line(atom: Atom, res: Residue, serial: int, record: Optional[str] = None) -> str:
    """ATOM/HETATM line (or, with ``record`` 'ANISOU', the ANISOU line of the atom)."""
    element = _element(atom, res.name)
    body = atom.tail.ljust(49)[:49] if record is None else atom.anisou.ljust(53)[:49]
    return (f"{record or atom.record:<6}{serial % 100000:>5} {_atom_name_field(atom.name, element)}"
            f"{(atom.alt or ' ')[:1]}{_id_columns(res)}{body}{element:>2}{atom.charge:<2}")


def _seqres_insert_index(header: List[str]) -> int:
    for k, line in enumerate(header):
        if line[:6].strip().upper().startswith(_SEQRES_ORDER):
            return k
    return len(header)


def write_pdb(structure: Structure, strip_conect: bool = True, seqres: bool = False) -> str:
    """Format a structure: serials renumbered (TER takes one), TER after every segment, END.

    ``strip_conect``: CONECT records should be stripped until final deposition processing; when kept they
    are remapped to the new serials.  ``seqres``: (re)generate the SEQRES records in the header.
    """
    header = list(structure.header)
    if seqres:
        header = [line for line in header if not line.startswith("SEQRES")]
        k = _seqres_insert_index(header)
        header[k:k] = seqres_records(structure)
    out = list(header)
    if structure.model:
        out.append(structure.model)
    serial = 0
    serial_map: Dict[int, List[int]] = {}
    for seg in structure.segments:
        if not seg:
            continue
        for res in seg:
            for atom in res.atoms:
                serial += 1
                if atom.serial is not None:
                    serial_map.setdefault(atom.serial, []).append(serial)
                out.append(_atom_line(atom, res, serial))
                if atom.anisou:
                    out.append(_atom_line(atom, res, serial, "ANISOU"))
        serial += 1
        out.append(f"TER   {serial % 100000:>5}      {_id_columns(seg[-1])}")
    if structure.model:
        out.append("ENDMDL")
    out.extend(structure.trailer)
    if not strip_conect:
        out.extend(_remap_conect(structure.conect, serial_map))
    out.append("END")
    return "\n".join(line.rstrip() for line in out) + "\n"


def _remap_conect(lines: Sequence[str], serial_map: Dict[int, List[int]]) -> List[str]:
    if any(len(v) > 1 for v in serial_map.values()):
        return []  # ambiguous input serials (wrapped numbering): CONECT cannot be remapped
    out = []
    for line in lines:
        fields = [line.ljust(31)[i:i + 5] for i in range(6, 31, 5)]
        try:
            ids = [int(f) for f in fields if f.strip()]
        except ValueError:
            continue
        if not ids or ids[0] not in serial_map:
            continue
        new = [serial_map[ids[0]][0]] + [serial_map[i][0] for i in ids[1:] if i in serial_map]
        if len(new) > 1:
            out.append("CONECT" + "".join(f"{i % 100000:>5}" for i in new))
    return out


# --- chirality ----------------------------------------------------------------------------------------

def chiral_volume(centre, a, b, c) -> float:
    """Signed volume a.(b x c) of the vectors from ``centre`` to atoms 1-3 (CIF _chem_comp_chir sign)."""
    va, vb, vc = ([float(p[i]) - float(centre[i]) for i in range(3)] for p in (a, b, c))
    cross = (vb[1] * vc[2] - vb[2] * vc[1], vb[2] * vc[0] - vb[0] * vc[2], vb[0] * vc[1] - vb[1] * vc[0])
    return va[0] * cross[0] + va[1] * cross[1] + va[2] * cross[2]


def sugar_chirality(coords: Dict[str, Vec], base: Optional[str] = None
                    ) -> Tuple[Optional[str], Dict[str, float]]:
    """Chirality call "D", "L" or None (no usable centre, or centres that disagree) from the C4', C3', C1'
    volumes; a nearly planar centre (below ``MIN_CHIRAL_VOLUME``) also prevents a call.

    ``coords``: atom name -> xyz of one residue; ``base`` selects the glycosidic nitrogen (N9 for A/G,
    N1 for C/U/T; unknown: N9 if present, else N1).  Returns the call and the volumes by centre.
    """
    if base in ("A", "G") or (base not in ("C", "U", "T") and "N9" in coords):
        glyco = "N9"
    else:
        glyco = "N1"
    volumes: Dict[str, float] = {}
    votes = Counter()
    for centre, partners, l_sign in L_CHIRAL_SIGNS:
        names = [glyco if p == "N" else p for p in partners]
        if centre not in coords or any(n not in coords for n in names):
            continue
        v = chiral_volume(coords[centre], *(coords[n] for n in names))
        volumes[centre] = v
        if abs(v) >= MIN_CHIRAL_VOLUME:
            votes["L" if v * l_sign > 0 else "D"] += 1
        else:
            votes["flat"] += 1  # a nearly planar centre: distorted geometry
    if not (votes["L"] or votes["D"]) or (votes["L"] and votes["D"]) or votes["flat"]:
        return None, volumes  # nothing usable, centres that disagree, or a flat centre: no call
    return ("L" if votes["L"] > votes["D"] else "D"), volumes


Range = Tuple[str, Optional[int], Optional[int]]


def parse_ranges(spec) -> List[Range]:
    """Residue ranges -> [(chain, start, end)].

    Accepts ``"C:45-73, A:1-5"``, ``"C:45"`` (one residue), ``"C"`` (whole chain), ``":45-73"`` (blank
    chain), a list of such strings, or (chain, start, end) tuples.
    """
    if spec is None:
        return []
    items = re.split(r"[,;\s]+", spec.strip()) if isinstance(spec, str) else list(spec)
    out: List[Range] = []
    for item in items:
        if isinstance(item, (tuple, list)):
            chain, a, b = (list(item) + [None, None])[:3]
            chain = str(chain).strip()
        elif not item:
            continue
        else:
            m = re.match(r"^([^:\s]?)(?::(-?\d+)(?:-(-?\d+))?)?$", item)
            if not m:
                raise PDBError(f"cannot read residue range {item!r} (expected e.g. C:45-73, C:45 or C)")
            chain, a, b = m.groups()
        start = None if a is None else int(a)
        end = start if b is None else int(b)
        out.append((chain, start, end))
    return out


def _in_ranges(res: Residue, ranges: Sequence[Range]) -> bool:
    for chain, start, end in ranges:
        if res.chain == chain and (start is None or start <= res.resseq <= end):
            return True
    return False


def _listing(items: Sequence[str], limit: int = 8) -> str:
    shown = ", ".join(items[:limit])
    return shown + (f" and {len(items) - limit} more" if len(items) > limit else "")


def _residue_ranges(residues: Sequence[Residue]) -> str:
    """'C:45-57, C:64-73' style summary of residues (consecutive numbers merged per chain)."""
    parts: List[List] = []
    for r in residues:
        if parts and parts[-1][0] == r.chain and r.resseq == parts[-1][2] + 1:
            parts[-1][2] = r.resseq
        else:
            parts.append([r.chain, r.resseq, r.resseq])
    return ", ".join(f"{c or '_'}:{a}" + (f"-{b}" if b != a else "") for c, a, b in parts)


def _volumes_text(vols: Dict[str, float]) -> str:
    return ", ".join(f"{k} {v:+.2f}" for k, v in vols.items())


def assign_chirality(structure: Structure, detect: bool = True, l_ranges=None) -> List[str]:
    """Set ``Residue.chirality`` of every nucleotide: override ranges (forced L) > detection > name."""
    ranges = parse_ranges(l_ranges)
    notes: List[str] = []
    l_in_d_names, d_in_l_names, undetermined, forced_d, l_dna = [], [], [], [], []
    counts = Counter()
    in_range = 0
    for res in structure.residues():
        nt = _nucleotide(res)
        if nt is None:
            continue
        call, vols = sugar_chirality(res.coords(), nt.base) if detect else (None, {})
        forced = bool(ranges) and _in_ranges(res, ranges)
        in_range += forced
        if forced:
            res.chirality = "L"
            if call == "D":
                forced_d.append(f"{res.label} ({_volumes_text(vols)})")
        else:
            res.chirality = call or nt.chirality
            if detect and call is None:
                undetermined.append(res.label)
            elif call == "L" and nt.chirality == "D":
                l_in_d_names.append(res)
            elif call == "D" and nt.chirality == "L":
                d_in_l_names.append(f"{res.label} ({_volumes_text(vols)})")
        if nt.sugar == "dna" and res.chirality == "L":
            l_dna.append(res.label)
        counts[res.chirality] += 1
    source = "sugar geometry" if detect else "residue names"
    notes.append(f"Chirality ({source}): {counts['L']} L and {counts['D']} D nucleotides"
                 + (f", {in_range} forced L by ranges" if ranges else ""))
    if l_in_d_names:
        notes.append(f"L sugars in {len(l_in_d_names)} residues with D names (e.g. cpptraj output), L names "
                     f"when a scheme is applied: {_residue_ranges(l_in_d_names)}")
    if d_in_l_names:
        notes.append("WARNING: L-named residues with D sugar geometry (D names when a scheme is applied; "
                     "force L with ranges if intended): " + _listing(d_in_l_names))
    if forced_d:
        notes.append("WARNING: forced L although the sugar geometry is D: " + _listing(forced_d))
    if undetermined:
        notes.append("Chirality undetermined (missing sugar atoms or centres that disagree), using the residue "
                     "name for "
                     + _listing(undetermined))
    if l_dna:
        notes.append("L-DNA residues have no name in the supported schemes and keep their names: "
                     + _listing(l_dna))
    if ranges and not in_range:
        notes.append("WARNING: the L ranges match no nucleotide")
    return notes


# --- residue names ------------------------------------------------------------------------------------

def _terminals(seg: Sequence[Residue]) -> Dict[int, str]:
    """id(residue) -> terminal state over the nucleotides of a segment."""
    nts = [r for r in seg if _nucleotide(r)]
    if len(nts) == 1:
        return {id(nts[0]): "N"}
    out = {id(r): "" for r in nts}
    if nts:
        out[id(nts[0])], out[id(nts[-1])] = "5", "3"
    return out


def convert_names(structure: Structure, scheme: str, only=None, terminals: bool = True) -> List[str]:
    """Rename RNA residues to ``scheme``; ``only`` (ranges) limits the renaming to those residues.

    ``terminals=False`` writes chain ends with the internal names (LA rather than LA5/LA3, A rather than
    A5), for leap set-ups whose libraries map terminal residues with addPdbResMap."""
    if scheme not in SCHEMES:
        raise PDBError(f"unknown residue name scheme {scheme!r} (choose from {', '.join(SCHEMES)})")
    ranges = parse_ranges(only)
    changes = Counter()
    term_mismatch, phosphate, single_l, kept_phosphate = [], [], [], []
    for seg in structure.segments:
        terms = _terminals(seg)
        for res in seg:
            nt = _nucleotide(res)
            if nt is None or nt.sugar != "rna" or (ranges and not _in_ranges(res, ranges)):
                continue
            term = seg_term = terms[id(res)] if terminals else ""
            chir = res.chirality or nt.chirality
            if scheme == "ccd" and chir == "L" and term in ("5", "N") and any(a.name == "P" for a in res.atoms):
                # 0A5 ... stand for a 5'-OH end; a 5' end with its phosphate is the CCD component itself
                term = ""  # (a single residue with a 5'-phosphate and a 3'-OH is the component, too)
                kept_phosphate.append(res.label)
            new = scheme_name(nt.base, chir, term, scheme)
            single_l_unit = seg_term == "N" and chir == "L" and nt.terminal == "5"  # LA5 stands in for LAN
            if terminals and nt.terminal and nt.terminal != seg_term and not single_l_unit:
                term_mismatch.append(res.label)  # the old name shows the misplaced terminal
            if new != res.name:
                changes[(res.name, new)] += 1
                res.name = new
            if scheme == "amber" and terms[id(res)] in ("5", "N") and any(a.name == "P" for a in res.atoms):
                phosphate.append(res.label)
            if scheme == "amber" and term == "N" and chir == "L":
                single_l.append(res.label)
    notes = []
    if changes:
        summary = ", ".join(f"{a}->{b} x{n}" for (a, b), n in sorted(changes.items(), key=lambda kv: -kv[1]))
        notes.append(f"Renamed {sum(changes.values())} residues to the {scheme} scheme: {summary}")
    else:
        notes.append(f"Residue names already follow the {scheme} scheme")
    if term_mismatch:
        notes.append("Terminal names that are not at a segment end (missing TER?): "
                     + _listing(term_mismatch))
    if phosphate:
        notes.append("WARNING: 5'-terminal residues carry P/OP1/OP2, which the Amber 5' units do not have "
                     "(tleap cannot type them; remove them with --strip-5p or cap the end): " + _listing(phosphate))
    if kept_phosphate:
        notes.append("5'-terminal L residues with a 5'-phosphate keep the CCD component names (0A ..., which "
                     "include P/OP1/OP2/OP3; the 5' variants 0A5 ... have a 5'-OH): " + _listing(kept_phosphate))
    if single_l:
        notes.append("nucleic12_LRNA.lib has no unit for a single-residue L segment; wrote the 5' unit "
                     "(its O3' lacks a hydrogen): " + _listing(single_l))
    return notes


def to_lrna(source, ranges, scheme: str = "amber", strip_conect: bool = True,
            terminals: bool = True) -> "CleanResult":
    """to_lrna.sh replacement: force L for ``ranges`` and rename only those residues (terminals kept unless
    ``terminals`` is False, which reproduces to_lrna.sh's plain LG/LA/LC/LU)."""
    s = read_pdb(source)
    notes = list(s.notes)
    notes += assign_chirality(s, detect=True, l_ranges=ranges)
    notes += convert_names(s, scheme, only=ranges, terminals=terminals)
    return CleanResult(write_pdb(s, strip_conect=strip_conect), notes, None, s)


def set_records(structure: Structure, policy: str) -> List[str]:
    """Apply a record policy: 'keep', 'rosetta' (polymer residues as ATOM) or 'pdb' (non-standard HETATM)."""
    if policy not in RECORD_POLICIES:
        raise PDBError(f"unknown record policy {policy!r} (choose from {', '.join(RECORD_POLICIES)})")
    if policy == "keep":
        return []
    changed = Counter()
    for res in structure.residues():
        if policy == "rosetta":
            record = "ATOM" if _is_polymer(res) else None
        else:
            record = "ATOM" if res.name in _STANDARD_NAMES else "HETATM"
        if record and any(a.record != record for a in res.atoms):
            changed[(res.name, record)] += 1
            for a in res.atoms:
                a.record = record
    if not changed:
        return []
    summary = ", ".join(f"{name}->{rec} x{n}" for (name, rec), n in sorted(changed.items()))
    return [f"Record policy {policy}: {summary}"]


# --- chains -------------------------------------------------------------------------------------------

def _resegment(structure: Structure) -> None:
    """Split segments wherever the chain id changes (a TER is written between chains)."""
    segments = []
    for seg in structure.segments:
        cur: List[Residue] = []
        for res in seg:
            if cur and cur[-1].chain != res.chain:
                segments.append(cur)
                cur = []
            cur.append(res)
        if cur:
            segments.append(cur)
    structure.segments = segments


def _polymer_segments(structure: Structure) -> List[List[Residue]]:
    return [seg for seg in structure.segments if not all(_is_solvent(r) for r in seg)]


def restore_chains(structure: Structure, reference=None, chain_ids: Optional[Sequence[str]] = None,
                   copy_numbers: bool = False) -> List[str]:
    """Restore chain ids (cpptraj output has none).

    ``reference``: copy ids residue by residue from a reference PDB (path, text or Structure); residues
    (water and ions excluded) are matched in order and compared by base letter, ignoring D/L and terminal
    naming; ``copy_numbers`` also copies residue numbers and insertion codes.  Otherwise ``chain_ids``
    gives one id per TER segment, else segments get A, B, C ...
    """
    if reference is not None:
        ref = read_pdb(reference)
        mine = [r for r in structure.residues() if not _is_solvent(r)]
        theirs = [r for r in ref.residues() if not _is_solvent(r)]
        if len(mine) != len(theirs):
            raise PDBError(f"chain restore: the structure has {len(mine)} residues but the reference has "
                           f"{len(theirs)} (water and ions excluded)")
        bad = [f"#{k + 1}: {a.label} vs reference {b.label}" for k, (a, b) in enumerate(zip(mine, theirs))
               if _code(a) != _code(b) or (_code(a) == "X" and a.name.upper() != b.name.upper())]
        if bad:
            raise PDBError(f"chain restore: sequence differs from the reference at {len(bad)} residue(s): "
                           + _listing(bad, 5))
        for a, b in zip(mine, theirs):
            a.chain = b.chain
            if copy_numbers:
                a.resseq, a.icode = b.resseq, b.icode
        _resegment(structure)
        return [f"Restored chain ids{' and residue numbers' if copy_numbers else ''} of {len(mine)} residues "
                f"from the reference (chains {', '.join(c or '_' for c in ref.chain_ids())})"]
    segments = _polymer_segments(structure)
    if chain_ids is not None:
        ids = [c.strip() for c in chain_ids]
        if len(ids) != len(segments):
            raise PDBError(f"chain restore: {len(ids)} chain ids given for {len(segments)} TER segments")
        if any(len(c) > 1 for c in ids):
            raise PDBError("chain ids must be single characters")
    else:
        if len(segments) > len(CHAIN_IDS):
            raise PDBError(f"{len(segments)} segments but only {len(CHAIN_IDS)} single-character chain ids")
        ids = list(CHAIN_IDS[:len(segments)])
    for seg, c in zip(segments, ids):
        for res in seg:
            res.chain = c
    _resegment(structure)
    return [f"Assigned chain ids {', '.join(c or '_' for c in ids)} to {len(segments)} TER segments"]


def _free_ids(structure: Structure) -> Iterator[str]:
    used = set(structure.chain_ids())
    for c in CHAIN_IDS:
        if c not in used:
            used.add(c)
            yield c
    raise PDBError("no unused single-character chain id left")


def split_chains(structure: Structure, at=None) -> List[str]:
    """Give new unused chain ids to every further TER segment of a chain (``at`` None), or split chains
    at residues (``at``: e.g. ``"C:45"`` - residue 45 and the rest of chain C get a new id and a TER)."""
    free = _free_ids(structure)
    notes = []
    if at is None:
        seen = set()
        for seg in _polymer_segments(structure):
            chain = seg[0].chain
            if chain in seen:
                new = next(free)
                notes.append(f"Segment starting at {seg[0].label} -> chain {new}")
                for res in seg:
                    res.chain = new
            seen.add(chain)
        return notes or ["No chain spans several TER segments; nothing split"]
    order = {id(r): n for n, r in enumerate(structure.residues())}
    targets = []  # residues to split before, found by their input chain and number, in file order
    for chain, start, _ in parse_ranges(at):
        if start is None:
            raise PDBError(f"split position needs a residue number, e.g. {chain or ''}:45")
        res = next((r for r in structure.residues() if r.chain == chain and r.resseq == start), None)
        if res is None:
            raise PDBError(f"split position {chain}:{start} not found")
        targets.append(res)
    for target in sorted(targets, key=lambda r: order[id(r)]):
        chain, start = target.chain, target.resseq  # the chain may already be a new id from an earlier split
        i, k = next((i, k) for i, seg in enumerate(structure.segments) for k, r in enumerate(seg) if r is target)
        seg = structure.segments[i]
        if k > 0:
            structure.segments[i:i + 1] = [seg[:k], seg[k:]]
            i += 1
        new = next(free)
        for later in structure.segments[i:]:
            for res in later:
                if res.chain == chain:
                    res.chain = new
        notes.append(f"Split chain {chain or '_'} at residue {start}: the rest is chain {new}")
    return notes


def combine_chains(structure: Structure, chains: Sequence[str], into: Optional[str] = None) -> List[str]:
    """Relabel ``chains`` as one chain (``into`` or the first listed); TER segments are kept."""
    chains = [c.strip() for c in chains]
    if not chains:
        return []
    target = chains[0] if into is None else into.strip()
    missing = [c for c in chains if c not in structure.chain_ids()]
    n = 0
    for res in structure.residues():
        if res.chain in chains:
            res.chain = target
            n += 1
    listed = ", ".join(c or "_" for c in chains)
    notes = [f"Combined chains {listed} into chain {target or '_'} ({n} residues)"]
    if missing:
        notes.append("WARNING: chains not present: " + ", ".join(missing))
    keys = Counter((r.chain, r.resseq, r.icode) for r in structure.residues() if r.chain == target)
    dup = sorted(k[1] for k, v in keys.items() if v > 1)
    if dup:
        notes.append(f"WARNING: chain {target or '_'} now has duplicate residue numbers ("
                     + _listing([str(d) for d in dup]) + "); renumber the residues")
    return notes


def order_chains(structure: Structure, order: Optional[Sequence[str]] = None) -> List[str]:
    """Reorder segments by chain: listed chains first in ``order`` (alphabetical when None); water and
    ion segments stay at the end in their original order."""
    before = structure.chain_ids()
    rank = {c: k for k, c in enumerate(order)} if order is not None else None

    def key(item):
        k, seg = item
        if all(_is_solvent(r) for r in seg):
            return (2, 0, "", k)
        c = seg[0].chain
        if rank is None:
            return (0, c == "", c, k)
        return (0, 0, "", rank[c]) if c in rank else (1, 0, "", k)

    structure.segments = [seg for _, seg in sorted(enumerate(structure.segments), key=key)]
    after = structure.chain_ids()
    return [] if after == before else [f"Chain order: {' '.join(c or '_' for c in after)}"]


# --- numbering and sequence ---------------------------------------------------------------------------

PHOSPHATE_5P = ("P", "OP1", "OP2", "OP3", "O1P", "O2P", "O3P")


def strip_5prime_phosphates(structure: Structure) -> List[str]:
    """Remove P/OP1/OP2/OP3 from the first nucleotide of every TER segment (Amber's 5'-terminal units
    have a 5'-OH).  A residue that has hydrogens gets HO5' where the phosphorus was (0.96 Å from O5')."""
    changed, skipped, emptied = [], [], []
    for seg in structure.segments:
        res = next((r for r in seg if not _is_solvent(r)), None)
        while res is not None and res.atoms and all(a.name.replace("*", "'") in PHOSPHATE_5P for a in res.atoms):
            emptied.append(res.label)  # a lone 5'-phosphate modelled as its own residue
            seg.remove(res)
            res = next((r for r in seg if not _is_solvent(r)), None)
        if res is None:
            continue
        if _nucleotide(res) is None:
            skipped.append(res.label)  # the segment starts with something else: its phosphates stay
            continue
        p_atoms = [a for a in res.atoms if a.name.replace("*", "'") in PHOSPHATE_5P]
        if not p_atoms:
            continue
        p = next((a for a in p_atoms if a.name == "P"), None)
        o5 = next((a for a in res.atoms if a.name.replace("*", "'") == "O5'"), None)
        has_h = any(a.name.lstrip("0123456789").startswith("H") for a in res.atoms)
        res.atoms = [a for a in res.atoms if a not in p_atoms]
        if has_h and p is not None and o5 is not None and not any(a.name == "HO5'" for a in res.atoms):
            o, d = o5.xyz, [p.xyz[k] - o5.xyz[k] for k in range(3)]
            n = math.sqrt(sum(v * v for v in d)) or 1.0
            h = [o[k] + 0.96 * d[k] / n for k in range(3)]
            tail = f"   {h[0]:8.3f}{h[1]:8.3f}{h[2]:8.3f}" + o5.tail[27:]
            res.atoms.insert(res.atoms.index(o5) + 1, Atom(o5.record, "HO5'", "", tail, "H"))
        changed.append(res.label)
    notes = []
    if changed:
        notes.append(f"Removed the 5'-terminal phosphate of {len(changed)} residue(s): " + _listing(changed))
    if emptied:
        notes.append(f"Removed {len(emptied)} 5'-terminal residue(s) that held only a phosphate: " + _listing(emptied))
    if skipped:
        notes.append("5' phosphates kept where a segment does not start with a recognised nucleotide: "
                     + _listing(skipped))
    return notes


def renumber_residues(structure: Structure, start: int = 1, per_chain: bool = False,
                      keep_gaps: bool = False) -> List[str]:
    """Renumber residues in file order, continuously across chains (each chain continues after the last
    number) or restarting at ``start`` per chain; ``keep_gaps`` keeps the numbering gaps (and insertion
    codes) inside a chain, otherwise numbers are sequential and insertion codes are cleared."""
    last: Dict[Optional[str], int] = {}
    prev: Optional[Residue] = None
    prev_old: Tuple[int, str] = (0, "")
    cleared = disorder = 0
    for res in structure.residues():
        key = res.chain if per_chain else None
        old = (res.resseq, res.icode)
        same_chain = prev is not None and prev.chain == res.chain
        if key not in last:
            new = start
        elif keep_gaps and same_chain and old > prev_old:
            new = last[key] + max(0, old[0] - prev_old[0])
        else:
            new = last[key] + 1
            disorder += keep_gaps and same_chain
        if not keep_gaps and res.icode:
            cleared += 1
            res.icode = ""
        prev, prev_old = res, old
        res.resseq = last[key] = new
    mode = "per chain" if per_chain else "continuously"
    notes = [f"Renumbered residues {mode} from {start}" + (" (gaps kept)" if keep_gaps else "")]
    if cleared:
        notes.append(f"Cleared {cleared} insertion codes")
    if disorder:
        notes.append(f"WARNING: {disorder} residues out of order; numbered sequentially there")
    if max((r.resseq for r in structure.residues()), default=0) > 9999:
        notes.append("WARNING: residue numbers above 9999 wrap in the PDB columns")
    return notes


def _sequence_chains(structure: Structure) -> Dict[str, List[Residue]]:
    """Chain -> residues of the polymer segments (water and ions excluded), in order."""
    chains: Dict[str, List[Residue]] = {}
    for seg in _polymer_segments(structure):
        if not any(_is_polymer(r) for r in seg):
            continue  # ligand-only segment
        for res in seg:
            if not _is_solvent(res):
                chains.setdefault(res.chain, []).append(res)
    return chains


def seqres_records(structure: Structure) -> List[str]:
    """SEQRES records (13 residue names per line) of the modelled residues, names as written."""
    out = []
    for chain, residues in _sequence_chains(structure).items():
        names = [r.name for r in residues]
        for k in range(0, len(names), 13):
            body = " ".join(f"{n:>3}" for n in names[k:k + 13])
            out.append(f"SEQRES {k // 13 + 1:>3} {(chain or ' ')[:1]} {len(names):>4}  {body}")
    return out


def _compress(numbers: Sequence[int]) -> str:
    parts = []
    for n in numbers:
        if parts and n == parts[-1][1] + 1:
            parts[-1][1] = n
        else:
            parts.append([n, n])
    return ", ".join(str(a) if a == b else f"{a}-{b}" for a, b in parts)


def fasta(structure: Structure, title: str = "structure") -> str:
    """FASTA text, one entry per chain; L-nucleotides are noted in the header line."""
    out = []
    for chain, residues in _sequence_chains(structure).items():
        nts = [r for r in residues if _nucleotide(r)]
        sugars = {_nucleotide(r).sugar for r in nts}
        n_amino = sum(r.name.upper() in _AMINO or r.name.upper() in _AMINO_VARIANTS for r in residues)
        if nts and not n_amino:
            kind = "RNA" if sugars == {"rna"} else "DNA" if sugars == {"dna"} else "nucleic acid"
        else:
            kind = "protein" if n_amino and not nts else "polymer"
        l_res = [r for r in nts if _chirality(r) == "L"]
        if l_res and len(l_res) == len(nts):
            kind = "L-" + kind
        elif l_res:
            kind += f", L-nucleotides at {_compress([r.resseq for r in l_res])}"
        out.append(f">{title}_{chain or '_'} chain {chain or '_'} | {kind} | {len(residues)} residues")
        out.append("".join(_code(r) for r in residues))
    return "\n".join(out) + "\n"


# --- convenience entry point --------------------------------------------------------------------------

@dataclass
class CleanResult:
    text: str  # the PDB text
    notes: List[str]  # human-readable change notes
    fasta: Optional[str]  # FASTA text when requested
    structure: Structure


def clean(source, *, scheme: Optional[str] = None, records: Optional[str] = None, strip_conect: bool = True,
          detect_l: bool = True, l_ranges=None,
          chain_reference=None, chain_ids: Optional[Sequence[str]] = None, assign_chains: bool = False,
          copy_numbers: bool = False,
          split_segments: bool = False, split_at=None, combine: Optional[Sequence[str]] = None,
          combine_into: Optional[str] = None, sort_chains: Union[bool, Sequence[str]] = False,
          renumber: bool = False, start: int = 1, per_chain: bool = False, keep_gaps: bool = False,
          seqres: bool = False, make_fasta: bool = False, title: Optional[str] = None,
          terminal_names: bool = True, strip_5p: bool = False) -> CleanResult:
    """Read, fix and write a PDB file in one call; returns the text plus change notes.

    Order: parse (ChimeraX repair, first model) -> chain restore (``chain_reference`` / ``chain_ids`` /
    ``assign_chains``) -> chirality (detection unless ``detect_l`` is False; ``l_ranges`` force L and use
    the input residue numbers with the restored chain ids) -> split / combine / ``sort_chains`` (True =
    alphabetical, or an explicit order) -> 5'-phosphate removal (``strip_5p``) -> residue names (``scheme``) ->
    renumbering -> record policy
    (``records``; None follows the scheme: 'rosetta' for rosetta, 'pdb' for pdb, else 'keep') -> writing.
    """
    s = read_pdb(source)
    notes = list(s.notes)
    stale = [line for line in s.header if line.startswith("REMARK 950")]
    if stale:  # base-pair annotation of the old names / numbers (mdatk.annotate)
        s.header = [line for line in s.header if not line.startswith("REMARK 950")]
        notes.append(f"Dropped {len(stale)} REMARK 950 base-pair annotation lines (they name the input "
                     "residues); annotate the cleaned file again if needed")
    if chain_reference is not None or chain_ids is not None or assign_chains:
        notes += restore_chains(s, chain_reference, chain_ids, copy_numbers)
    if detect_l or l_ranges:
        notes += assign_chirality(s, detect=detect_l, l_ranges=l_ranges)
    if split_segments:
        notes += split_chains(s)
    if split_at:
        notes += split_chains(s, split_at)
    if combine:
        notes += combine_chains(s, combine, combine_into)
    if sort_chains:
        notes += order_chains(s, None if sort_chains is True else list(sort_chains))
    if strip_5p:
        notes += strip_5prime_phosphates(s)
    if scheme:
        notes += convert_names(s, scheme, terminals=terminal_names)
    if renumber:
        notes += renumber_residues(s, start, per_chain, keep_gaps)
    policy = records or {"rosetta": "rosetta", "pdb": "pdb"}.get(scheme or "", "keep")
    notes += set_records(s, policy)
    if s.conect:
        notes.append(f"Stripped {len(s.conect)} CONECT records" if strip_conect else
                     "Kept CONECT records (remapped to the new atom serials)")
    n_atoms = sum(len(r.atoms) for r in s.residues())
    if n_atoms + len(s.segments) > 99999:
        notes.append("WARNING: more than 99999 atoms; serial numbers wrap in the PDB columns")
    if seqres:
        notes.append(f"Wrote SEQRES records for {len(_sequence_chains(s))} chain(s)")
    name = title or "structure"
    is_path = isinstance(source, os.PathLike) or (isinstance(source, str) and "\n" not in source)
    if title is None and is_path:
        name = os.path.splitext(os.path.basename(str(source)))[0]
    text = write_pdb(s, strip_conect=strip_conect, seqres=seqres)
    return CleanResult(text, notes, fasta(s, name) if make_fasta else None, s)


# --- command line (``mdatk pdb`` and the copies embedded in generated scripts) -------------------------

def summary(structure: Structure) -> Dict[str, object]:
    """Chains, sequences, residue names and chirality calls of a parsed structure (for 'check')."""
    chains: Dict[str, Dict[str, object]] = {}
    for chain, residues in _sequence_chains(structure).items():
        calls = Counter(_chirality(r) or "?" for r in residues if _nucleotide(r))
        chains[chain or "_"] = {
            "residues": len(residues),
            "numbers": f"{residues[0].resseq}-{residues[-1].resseq}" if residues else "",
            "sequence": "".join(_code(r) for r in residues),
            "names": sorted({r.name for r in residues}),
            "chirality": dict(calls),
        }
    other = Counter(r.name for r in structure.residues() if not _is_polymer(r))
    return {"chains": chains, "segments": len(structure.segments), "other_residues": dict(other),
            "conect_records": len(structure.conect), "notes": list(structure.notes)}


def _clean_arguments(p) -> None:
    p.add_argument("input", help="PDB file")
    p.add_argument("--output", "-o", help="output PDB (default: print to stdout)")
    p.add_argument("--scheme", choices=SCHEMES, help="residue names: amber (LA/LA5... for the L-RNA force "
                   "field), ccd (0A... for Coot/Phenix), rosetta (0A..., all ATOM), pdb")
    p.add_argument("--records", choices=RECORD_POLICIES, help="ATOM/HETATM policy (default: follows the "
                   "scheme)")
    p.add_argument("--no-terminal-names", action="store_true", help="chain ends keep the internal names "
                   "(LA, not LA5/LA3), for libraries that map terminals with addPdbResMap")
    p.add_argument("--keep-conect", action="store_true", help="keep CONECT records (stripped by default "
                   "until final deposition processing)")
    p.add_argument("--no-detect", action="store_true", help="do not detect L-nucleotides from the sugar "
                   "chirality (C1'/C3'/C4'); use the residue names")
    p.add_argument("--l-ranges", help="force L for these residues, e.g. 'C:45-73' or 'C' (whole chain)")
    p.add_argument("--chain-reference", help="copy chain ids residue by residue from this PDB (e.g. the "
                   "structure a cpptraj PDB was built from)")
    p.add_argument("--chain-ids", help="chain ids per TER segment, e.g. A,B,C")
    p.add_argument("--assign-chains", action="store_true", help="assign A, B, C ... per TER segment")
    p.add_argument("--copy-numbers", action="store_true", help="with --chain-reference: also copy the "
                   "residue numbers")
    p.add_argument("--split-segments", action="store_true", help="give every TER segment its own new "
                   "chain id")
    p.add_argument("--split-at", help="split a chain before these residues, e.g. 'A:31'")
    p.add_argument("--combine", help="merge these chains, e.g. A,B")
    p.add_argument("--combine-into", help="chain id of the merged chain (default: the first)")
    p.add_argument("--sort-chains", action="store_true", help="write the chains in alphabetical order")
    p.add_argument("--renumber", action="store_true", help="renumber residues continuously across chains")
    p.add_argument("--start", type=int, default=1, help="first residue number (default 1)")
    p.add_argument("--per-chain", action="store_true", help="restart the numbering in every chain")
    p.add_argument("--keep-gaps", action="store_true", help="keep numbering gaps when renumbering")
    p.add_argument("--seqres", action="store_true", help="write SEQRES records")
    p.add_argument("--strip-5p", action="store_true", help="remove the phosphate (P/OP1/OP2/OP3) of every "
                   "chain's 5'-terminal nucleotide, as Amber's 5' units expect")
    p.add_argument("--fasta", help="also write the sequences to this FASTA file")


def clean_options(args) -> Dict[str, object]:
    """Keyword arguments of :func:`clean` from parsed command-line arguments."""
    def ids(text):
        return [c.strip() for c in text.split(",") if c.strip()] if text else None

    return dict(scheme=args.scheme, records=args.records, strip_conect=not args.keep_conect,
                detect_l=not args.no_detect, l_ranges=args.l_ranges or None,
                chain_reference=args.chain_reference, chain_ids=ids(args.chain_ids),
                assign_chains=args.assign_chains, copy_numbers=args.copy_numbers,
                split_segments=args.split_segments, split_at=args.split_at or None,
                combine=ids(args.combine), combine_into=args.combine_into, sort_chains=args.sort_chains,
                renumber=args.renumber, start=args.start, per_chain=args.per_chain, keep_gaps=args.keep_gaps,
                seqres=args.seqres, make_fasta=bool(args.fasta), terminal_names=not args.no_terminal_names,
                strip_5p=args.strip_5p)


def add_commands(sub) -> list:
    """Register 'clean', 'check' and 'to-lrna' on an argparse sub-parser collection (returns the parsers)."""
    c = sub.add_parser("clean", help="repair, rename, re-chain, renumber and rewrite a PDB file")
    _clean_arguments(c)
    k = sub.add_parser("check", help="report chains, sequences, residue names and sugar chirality")
    k.add_argument("input", help="PDB file")
    k.add_argument("--l-ranges", help="residues to treat as L, e.g. 'C:45-73'")
    t = sub.add_parser("to-lrna", help="rename residues in ranges to L-RNA names (replaces to_lrna.sh)")
    t.add_argument("input", help="PDB file")
    t.add_argument("ranges", help="e.g. 'C:45-73' or 'A:1-5,A:9-12'")
    t.add_argument("--scheme", choices=SCHEMES, default="amber")
    t.add_argument("--no-terminal-names", action="store_true", help="plain LG/LA/LC/LU at chain ends, like "
                   "to_lrna.sh")
    t.add_argument("--output", "-o", help="output PDB (default: stdout)")
    return [c, k, t]


def run(args, stdout=None, stderr=None) -> Tuple[int, Dict[str, object]]:
    """Execute a parsed 'clean' / 'check' / 'to-lrna' command; returns (exit status, result)."""
    import json

    stdout = stdout or sys.stdout
    stderr = stderr or sys.stderr
    if args.pdb_command == "check":
        s = read_pdb(args.input)
        notes = assign_chirality(s, detect=True, l_ranges=args.l_ranges or None)
        info = summary(s)
        info["notes"] = list(info["notes"]) + notes
        stdout.write(json.dumps(info, indent=2) + "\n")
        return 0, info
    if args.pdb_command == "to-lrna":
        res = to_lrna(args.input, args.ranges, scheme=args.scheme, terminals=not args.no_terminal_names)
    else:
        res = clean(args.input, **clean_options(args))
    if args.output:
        with open(args.output, "w", encoding="utf-8") as fh:
            fh.write(res.text)
    else:
        stdout.write(res.text)
    fasta_path = getattr(args, "fasta", None)
    if fasta_path and res.fasta:
        with open(fasta_path, "w", encoding="utf-8") as fh:
            fh.write(res.fasta)
    for note in res.notes:
        stderr.write(f"pdbtools: {note}\n")
    files = [f for f in (args.output, fasta_path) if f]
    return 0, {"notes": res.notes, "files": files}


def main(argv: Optional[Sequence[str]] = None) -> int:
    import argparse

    parser = argparse.ArgumentParser(prog=os.path.basename(sys.argv[0]) or "pdbtools",
                                     description="PDB hygiene for D- and L-RNA structures.")
    add_commands(parser.add_subparsers(dest="pdb_command", required=True, metavar="command"))
    args = parser.parse_args(argv)
    try:
        return run(args)[0]
    except (OSError, PDBError) as exc:
        sys.stderr.write(f"error: {exc}\n")
        return 1


# --- settings written by MDATK ---


if __name__ == "__main__":
    sys.exit(main())
