"""MDBridge's deliberate differences in MDATK v0.0.4.1's chemistry layer (docs/specs/m6-plan.md
section 5, CH-01 ... CH-15, owner decisions D7, D11 ... D13), as test-side patches of MDATK, so that
the chem generators (make_pdbtools_cases.py, make_annotate_cases.py, make_builder_cases.py,
make_dict_cases.py, make_density_fit_cases.py) and cli_parity.py record MDBridge's behaviour.
Where a later MDATK fixed the item, its design is followed and named (MDATK's own code, ported).

apply_mdatk() patches the library functions (idempotent); apply_cli() also the command line (the
child of cli_parity.py calls both). Each patch is listed with its id:

  CH-01  residue numbers below -999 do not fit columns 23-26: PDBError "residue number -1002 (A)
         does not fit the PDB columns (-999 to 9999)" when the record is written (pdbtools
         _id_columns, so `pdb clean`, `build helix`, `build mutate` alike); numbers above 9999 keep
         MDATK's wrap. Reading a line MDATK wrote that way ("A-1234" and the coordinates one column
         to the right) says so: "line 1: residue number -1234 does not fit the PDB columns (-999 to
         9999); the later columns are shifted" instead of "line 1: unreadable coordinates"
  CH-02  build_helix notes "WARNING: residue numbers above 9999 wrap in the PDB columns" (the
         renumbering's warning) when a strand runs past 9999; appended after the build notes
  CH-03  combine_chains(into=...) with more than one character: PDBError("chain ids must be single
         characters"), as restore_chains
  CH-04  clean(strip_conect=False) with repeated input serials: the note says "WARNING: dropped N
         CONECT records: the input's atom serials repeat, so they cannot be remapped" instead of
         "Kept CONECT records (remapped to the new atom serials)"
  CH-05  a residue with an RNA name and a 2'-deoxy sugar (C2' and no O2') is DNA (MDATK 0.1.0,
         pdbtools.py _nucleotide)
  CH-06  a chain end next to a residue the schemes do not know (PSU, a 5'-GTP) bonded to it through
         O3'-P (<= 2.0 A) is no end: the nucleotide keeps its internal name (MDATK 0.1.0,
         pdbtools.py _terminals / _bonded_on, without its LINK records)
  CH-07  SEQRES names have three characters: a longer nucleotide name as its CCD component (also
         the Amber-style L-DNA LDA ...), any other as UNK, with notes after "Wrote SEQRES records"
         (MDATK 0.0.17, pdbtools.py _seqres_names; its note texts)
  CH-08  a mutation without a chain id reads the whole number: 15G = blank chain, residue 15 (MDATK
         0.1.0, builder.py _MUTATION_RE; bases A C G U as in 0.0.4.1)
  CH-09  (D11) a mutated base gets the builder's ideal hydrogens (template rules, the same fit as its
         heavy atoms, after them in _H_RULES order) when the old base had hydrogens; the note
         "X:5: replaced 2 hydrogen(s) of the old base (H8 H2) with 3 of the new base (H8 H1 H21 ...)
         at ideal positions" instead of "dropped ...; the new base has heavy atoms only"
  CH-10  non-finite resolutions and negative grid spacings are refused before any work: build fit
         (a NaN / infinite --resolution: "the resolution must be a positive number of Å"; zero and
         negative ones keep MDATK's "enter the map resolution (Å)"), density fit (NaN / infinite;
         negative ones keep the project check "Density map resolution must be positive."), density
         simulate (any resolution that is not a positive finite number; a negative or non-finite
         --spacing: "the grid spacing must be positive"; --spacing 0 keeps the default)
  CH-11  non-finite map voxels carry no density in a fit (RN-30's rule): the fit runs on a copy with
         them as 0, the density centre of a local start uses the finite voxels (RN-33's level), and
         the fit message ends "; N non-finite map voxel(s) read as no density"
  CH-12  density fit --save writes the project through the edit model (ED-02): the file's own paths as
         it had them (relative stays relative), only the density settings changed; a --map given
         relative is stored relative to the project folder; --save OTHER.mdatk.toml rebases the
         relative paths onto the other folder (Save As). Without --project MDATK's save is kept
  CH-13  the command line names a missing input file: "error: [Errno 2] No such file or directory:
         'x.pdb'" (a folder: "[Errno 21] Is a directory: 'x'") for pdb clean / check / to-lrna /
         annotate, build mutate / fit and --chain-reference; text given to the API is unchanged
  CH-14  (D9) annotate: 'auto' is the builtin engine; 'barnaba' runs the builtin engine with the note
         "WARNING: Barnaba is not available in MDBridge; using the builtin detector" (before the
         engine's notes); the stacking note without the pip hint
  CH-15  a reversed residue range (A:4-1) is an error: PDBError("residue range A:4-1 ends before it
         starts")
  CHD-01 (owner-approved) guanine's amino hydrogens as the wwPDB CCD names them: H21 trans to N1
         (N3 side), H22 cis to N1 (the Watson-Crick edge): the builder's _H_RULES["G"] (built and
         mutated G), the Watson-Crick-edge hydrogen of N2 in _PAIR_HBONDS (wc and reverse-wc G·C,
         reverse-Hoogsteen G·G: H22) and the 0G dictionary rule (COMPONENTS 0G / 0G5 / 0G3 hpos, from
         tools/gen_chem_tables.hydrogen_rules()); 0.1.1 not changed
  CHD-02 (owner-approved) the 4JRD nucleotide in the frame of the fitted 4JRD screw, with its twist and
         rise (_JRD_NUCLEOTIDE, JRD_TWIST, JRD_RISE from tools/gen_chem_tables.jrd_nucleotide(), the
         tables of src/chem/helix_tables.inc): every non-Watson-Crick build, and REMARK 951's hh TWIST
         41.7 (MDATK: 41.8); 0.1.1 not changed
  CH-16  (owner decision D17) the other pairings' strand closure (builder._close_levels) by
         continuation: from the uniform helix (every partner nucleoside on its strand's reference
         pose, closed by _helix_fit's torsions) the partner poses move to their _level_pose poses
         along the rigid transform (rotation vector about the reference C1' and C1' shift, scaled by
         t) in CLOSURE_STEPS steps; a step whose solve stops on its evaluation limit, or whose
         contact pairs need more than CLOSURE_ROUNDS_OK refreshes, is retried as two half steps
         (recursively, down to 1/320 of the interval). Every solve: MDATK's residuals at its last
         pass's contact weight 0.01, TRF with the exact subproblem (2-point differences grouped by
         MDATK's level-span sparsity, then dense), ftol = xtol = gtol = 1e-10; the contact pairs
         within 4.5 Å (MDATK: 3.8) refreshed before every solve, kept once scored, in (p, q) order,
         the solve repeated while the set grows; a solve's variables are the move from its start (TRF's
         first trust radius 1, not the start's norm). _pair_fit, _level_pose and _helix_fit converged
         (ftol = xtol = gtol = 1e-14). MDATK's four truncated passes ended in machine-dependent
         minima (up to 1.6 Å RMSD between kernels); 0.1.1 has the same closure. No SciPy private API:
         the column groups and the differences are written out (column_groups, fd_jacobian; SciPy's
         group_columns / approx_derivative give the same numbers)
  CH-17  the link torsions' term of the other pairings' fits is continuous (builder._link_residuals:
         MDATK's (1 + cos 3 theta) * 1.5 + wrap(theta - ref) / 120 jumps by 3 at theta = ref + 180,
         where its reference helices sit and its closures creep): the tilt TILT_A sin d +
         TILT_B sin 2d (d = theta - ref) has MDATK's slope at the reference and its values +-1 at the
         staggered positions +-120 degrees; the closure always uses it, the reference helices' fit
         as REFERENCE_TERM says (the reference geometry is a separable input; which one is right is
         open); 0.1.1 has MDATK's term
  D7     dict pubchem asks MDBRIDGE_PUBCHEM_URL (the PUG REST root, default MDATK's PUBCHEM_REST) when
         that variable is set: the local stand-in server of cli_parity.py
  D9(M4) retired: `render` / `preview` fit a shown map without a placement again when auto_fit is on
         (MDATK's own ensure_density_fit; improvements_render.py no longer patches it)
  D12    unknown residues pass through (v0.0.4.1) plus CH-06: nothing to patch
  D13    REMARK 950 MDBridge BASE-PAIR ANNOTATION / REMARK 951 MDBridge IDEAL ...; every other byte as
         MDATK writes it
  names  messages that name the program: 'mdatk pdb clean' / 'mdatk density fit --save' as
         'mdbridge ...' (NAMING)

Not patched (kept as MDATK writes them, docs/specs/m6-plan.md section 5 "Not bugs, kept"): the
dictionaries' header comment "renamed by MDATK dictionaries.py" (milestone 7 decides the embedded
copies, plan 4.5).
"""

from __future__ import annotations

import copy
import errno
import math
import os
import re
from collections import Counter

NAMING = [("'mdatk pdb clean'", "'mdbridge pdb clean'"),
          ("'mdatk density fit --save'", "'mdbridge density fit --save'")]

BARNABA_NOTE = "WARNING: Barnaba is not available in MDBridge; using the builtin detector"
STACKING_NOTE = "The builtin engine finds canonical base pairs only (no stacking)"
PUBCHEM_URL_VARIABLE = "MDBRIDGE_PUBCHEM_URL"
O3_P_BONDED = 2.0  # Å (MDATK 0.1.0)
_MUTATION_RE = re.compile(r"^(?:([^:\s]?):|([^:\s\d\-]?))(-?\d+)([A-Za-z]?)(?:>|:|->)?([ACGUacgu])$")

_applied = False
_applied_cli = False


def rename(text: str) -> str:
    for old, new in NAMING:
        text = text.replace(old, new)
    return text


def _ccd_component(info) -> str:
    """The CCD component of a nucleotide name (no terminal variant)."""
    from mdatk import pdbtools
    if info.sugar == "dna":
        return ("0D" if info.chirality == "L" else "D") + info.base
    return pdbtools.scheme_name(info.base, info.chirality, "", "ccd")


def apply_mdatk():
    global _applied
    if _applied:
        return
    _applied = True
    _patch_names()
    _patch_reference_data()
    _patch_pdbtools()
    _patch_annotate()
    _patch_builder()
    _patch_density()
    _patch_torsions()
    _patch_closure()


# ---------------------------------------------------------------------------------------------
# reference data (CHD-01, CHD-02: docs/DIFFERENCES.md "Chemistry reference data")

def chem_tables():
    """tools/gen_chem_tables.py as a module (pure Python, reads data/sources only)."""
    import importlib.util
    import sys

    name = "mdbridge_gen_chem_tables"
    if name not in sys.modules:
        path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "tools", "gen_chem_tables.py")
        spec = importlib.util.spec_from_file_location(name, os.path.normpath(path))
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        sys.modules[name] = module
    return sys.modules[name]


# the pair hydrogen bonds whose N2 donor is guanine's Watson-Crick-edge hydrogen (MDATK: H21)
_N2_WC_EDGE_PAIRS = (("wc", "G", "C"), ("reverse-wc", "G", "C"), ("reverse-hoogsteen", "G", "G"))


def _patch_reference_data():
    import dataclasses
    from functools import _lru_cache_wrapper

    from mdatk import builder as B
    from mdatk import dictionaries as D

    tables = chem_tables()
    # CHD-01: builder rules (nh2: parallel to the ring bond refs[1] -> refs[0], so ("C2", "N1") is trans
    # to N1 and ("C2", "N3") cis to N1); the H21, H22 order is kept
    rules = dict((h, (x, rule, refs)) for h, x, rule, refs in B._H_RULES["G"])
    assert rules["H21"][2] == ("C2", "N3") and rules["H22"][2] == ("C2", "N1"), "MDATK's G rules changed"
    B._H_RULES["G"] = tuple((h, x, rule, {"H21": ("C2", "N1"), "H22": ("C2", "N3")}.get(h, refs))
                            for h, x, rule, refs in B._H_RULES["G"])
    for pairing, a, b in _N2_WC_EDGE_PAIRS:
        for key in {(pairing, a, b), (pairing, b, a)}:
            B._PAIR_HBONDS[key] = tuple((x, y, side, "H22" if (h == "H21" and "N2" in (x, y)) else h)
                                        for x, y, side, h in B._PAIR_HBONDS[key])
    # CHD-01: the 0G dictionary rule, from the CCD (gen_chem_tables.hydrogen_rules)
    ccd_rules = {r[0]: r for r in tables.hydrogen_rules()["0G"] if r[1] == "N2"}
    for code in ("0G", "0G5", "0G3"):
        comp = D.COMPONENTS[code]
        D.COMPONENTS[code] = dataclasses.replace(
            comp, hpos=tuple(ccd_rules.get(r[0], r) if r[1] == "N2" else r for r in comp.hpos))
    # CHD-02: the 4JRD nucleotide in the fitted screw's frame
    twist, rise, rows, _frame = tables.jrd_nucleotide()
    assert [r[0] for r in rows] == [r[0] for r in B._JRD_NUCLEOTIDE], "4JRD atom names differ"
    B._JRD_NUCLEOTIDE = tuple(tuple(r) for r in rows)
    B.JRD_TWIST, B.JRD_RISE = twist, rise
    for value in vars(B).values():
        if isinstance(value, _lru_cache_wrapper):
            value.cache_clear()


# ---------------------------------------------------------------------------------------------
# pdbtools

def _patch_names():
    """NA-01 in the chem layer: the L-DNA names, the CCD's own atom names of 0DA / 0DT and the backbone
    guard of the base-pair detector (improvements_nucleotides.py, its secstruct part only: the PDB
    tools of MDBridge read residue names with io::parseResidueName)."""
    import improvements_nucleotides as na
    from mdatk import pdbtools, secstruct

    if not hasattr(pdbtools, "_mdbridge_original_parse"):
        pdbtools._mdbridge_original_parse = pdbtools.parse_residue_name
    pdbtools.parse_residue_name = na.parse_residue_name
    if getattr(secstruct, "_mdbridge_names", False) or na._applied:
        return  # improvements_nucleotides.apply_mdatk() ran (it must run first when both are applied)
    secstruct._mdbridge_names = True
    secstruct.parse_residue_name = na.parse_residue_name
    init = secstruct.PairDetector.__init__

    def detector_init(self, resnames, atom_index, cutoff=3.5, wobble=True):
        index = [na.standard_index(r, idx) for r, idx in zip(resnames, atom_index)]
        names = [r if na.has_backbone_atoms(idx) else "" for r, idx in zip(resnames, index)]
        init(self, names, index, cutoff=cutoff, wobble=wobble)

    secstruct.PairDetector.__init__ = detector_init


def _patch_pdbtools():
    from mdatk import pdbtools as pt

    def _nucleotide(res):  # CH-05 (MDATK 0.1.0)
        info = pt.parse_residue_name(res.name)
        if info is None or not any("'" in a.name or "*" in a.name or a.name == "P" for a in res.atoms):
            return None
        if info.sugar == "rna":
            names = {a.name.replace("*", "'") for a in res.atoms}
            if "C2'" in names and "O2'" not in names:
                info = info._replace(sugar="dna")
        return info

    pt._nucleotide = _nucleotide

    def _bonded_on(before, after):  # CH-06 (MDATK 0.1.0)
        o3, p = before.coords().get("O3'"), after.coords().get("P")
        return o3 is not None and p is not None and math.dist(o3, p) <= O3_P_BONDED

    def _phosphate_only(res):
        return bool(res.atoms) and all(a.name.replace("*", "'") in pt.PHOSPHATE_5P for a in res.atoms)

    def _ends_state(five, three):
        return "N" if five and three else "5" if five else "3" if three else ""

    def _terminals(seg):
        nts = [r for r in seg if pt._nucleotide(r)]
        out = {id(r): "" for r in nts}
        if not nts:
            return out
        first, last = nts[0], nts[-1]
        k5, k3 = (next(k for k, r in enumerate(seg) if r is end) for end in (first, last))
        five = not (k5 > 0 and _bonded_on(seg[k5 - 1], first))
        three = not (k3 + 1 < len(seg) and not _phosphate_only(seg[k3 + 1]) and _bonded_on(last, seg[k3 + 1]))
        if len(nts) == 1:
            out[id(first)] = _ends_state(five, three)
        else:
            out[id(first)], out[id(last)] = _ends_state(five, False), _ends_state(False, three)
        return out

    pt._terminals = _terminals

    parse_ranges = pt.parse_ranges

    def parse_ranges_checked(spec):  # CH-15
        out = parse_ranges(spec)
        for chain, start, end in out:
            if start is not None and end is not None and end < start:
                raise pt.PDBError(f"residue range {chain}:{start}-{end} ends before it starts")
        return out

    pt.parse_ranges = parse_ranges_checked

    id_columns = pt._id_columns

    def _id_columns(res):  # CH-01
        if res.resseq < -999:
            raise pt.PDBError(f"residue number {res.resseq} ({res.chain or '_'}) does not fit the PDB columns "
                              "(-999 to 9999)")
        return id_columns(res)

    pt._id_columns = _id_columns

    parse_atom = pt._parse_atom

    def _parse_atom(line, lineno):  # CH-01 when reading: what MDATK wrote for such numbers
        try:
            return parse_atom(line, lineno)
        except pt.PDBError as exc:
            padded = line.ljust(80)
            if (str(exc) == f"line {lineno}: unreadable coordinates" and re.fullmatch(r"-[0-9]{4}", padded[22:27])
                    and pt._coords_ok(padded[28:77])):
                raise pt.PDBError(f"line {lineno}: residue number {int(padded[22:27])} does not fit the PDB columns "
                                  "(-999 to 9999); the later columns are shifted") from None
            raise

    pt._parse_atom = _parse_atom

    combine = pt.combine_chains

    def combine_chains(structure, chains, into=None):  # CH-03
        if into is not None and len(into.strip()) > 1:
            raise pt.PDBError("chain ids must be single characters")
        return combine(structure, chains, into)

    pt.combine_chains = combine_chains

    def _seqres_names(structure):  # CH-07 (MDATK 0.0.17)
        out, mapped, unknown = {}, Counter(), []
        for chain, residues in pt._sequence_chains(structure).items():
            names = []
            for res in residues:
                name = res.name.strip()
                nt = pt._nucleotide(res)
                amber_l_dna = (nt is not None and nt.sugar == "dna" and nt.chirality == "L"
                               and name.upper()[:2] == "LD")
                if len(name) > 3 or amber_l_dna:
                    info = pt.parse_residue_name(name)
                    if info is not None:
                        new = _ccd_component(info)
                        mapped[(name, new)] += 1
                    else:
                        new = "UNK"
                        unknown.append(res.label)
                    name = new
                names.append(name)
            out[chain] = names
        notes = []
        if mapped:
            notes.append("SEQRES names are the CCD's three-character ones: wrote the CCD components of MDATK's "
                         "Amber-style L-DNA names and of other four-character names "
                         + ", ".join(f"{a}->{b} x{n}" for (a, b), n in sorted(mapped.items())))
        if unknown:
            notes.append("WARNING: SEQRES names have three characters: wrote UNK for the residue name(s) of "
                         + pt._listing(unknown))
        return out, notes

    pt._mdbridge_seqres_names = _seqres_names

    def seqres_records(structure):
        out = []
        for chain, names in _seqres_names(structure)[0].items():
            for k in range(0, len(names), 13):
                body = " ".join(f"{n:>3}" for n in names[k:k + 13])
                out.append(f"SEQRES {k // 13 + 1:>3} {(chain or ' ')[:1]} {len(names):>4}  {body}")
        return out

    pt.seqres_records = seqres_records

    clean = pt.clean
    kept = "Kept CONECT records (remapped to the new atom serials)"

    def clean_improved(source, **kw):
        res = clean(source, **kw)
        notes = res.notes
        if kept in notes:  # CH-04
            serials = Counter(a.serial for seg in res.structure.segments for r in seg for a in r.atoms
                              if a.serial is not None)
            if any(n > 1 for n in serials.values()):
                notes[notes.index(kept)] = (f"WARNING: dropped {len(res.structure.conect)} CONECT records: the "
                                            "input's atom serials repeat, so they cannot be remapped")
        if kw.get("seqres"):  # CH-07
            at = next(k for k, n in enumerate(notes) if n.startswith("Wrote SEQRES records"))
            notes[at + 1:at + 1] = _seqres_names(res.structure)[1]
        return res

    pt.clean = clean_improved


# ---------------------------------------------------------------------------------------------
# annotate

def _patch_annotate():
    from mdatk import annotate as an
    from mdatk import pdbtools as pt

    an._nucleotide = pt._nucleotide
    an.parse_residue_name = pt.parse_residue_name
    original = an.annotate
    header_old = f"{an.REMARK} MDATK BASE-PAIR ANNOTATION"
    header_new = f"{an.REMARK} MDBridge BASE-PAIR ANNOTATION"  # D13

    def annotate(source, pairs=True, stacking=True, dotbracket=True, engine="auto", cutoff=3.5, wobble=True):
        if engine not in an.ENGINES:
            raise ValueError(f"unknown annotation engine {engine!r} (choose from {', '.join(an.ENGINES)})")
        res = original(source, pairs=pairs, stacking=stacking, dotbracket=dotbracket, engine="builtin",
                       cutoff=cutoff, wobble=wobble)
        notes = [rename(n) for n in res.notes]
        notes = [STACKING_NOTE if n.startswith(STACKING_NOTE) else n for n in notes]
        if engine == "barnaba":  # CH-14
            at = next(k for k, n in enumerate(notes)
                      if n == STACKING_NOTE or n.startswith("Left ") or n.startswith("Residue numbering gaps")
                      or re.match(r"Annotated \d+ base pair", n))
            notes.insert(at, BARNABA_NOTE)
        res.notes = notes
        res.records = [header_new if r == header_old else r for r in res.records]
        res.text = res.text.replace(header_old + "\n", header_new + "\n", 1)
        return res

    an.annotate = annotate


# ---------------------------------------------------------------------------------------------
# builder

def _patch_builder():
    import numpy as np
    from mdatk import builder as B
    from mdatk import pdbtools as pt

    def parse_mutations(spec):  # CH-08 (MDATK 0.1.0)
        items = re.split(r"[,;\s]+", spec.strip()) if isinstance(spec, str) else list(spec)
        out = []
        for item in items:
            if not item:
                continue
            m = _MUTATION_RE.match(item.strip())
            if not m:
                raise B.BuildError(f"cannot read mutation {item!r} (expected e.g. A:5G - chain, residue number, "
                                   "new base)")
            with_colon, bare, num, icode, base = m.groups()
            mut = B.Mutation(with_colon or bare or "", int(num), icode, base.upper())
            if any((o.chain, o.resseq, o.icode) == (mut.chain, mut.resseq, mut.icode) for o in out):
                raise B.BuildError(f"residue {mut.where} is mutated twice")
            out.append(mut)
        if not out:
            raise B.BuildError("no mutations given (e.g. A:5G, B:20C)")
        return out

    B.parse_mutations = parse_mutations

    def _swap_base(res, base):  # MDATK 0.0.4.1's, with CH-09 (D11)
        from mdatk.spatial import kabsch

        info = pt.parse_residue_name(res.name)
        old = info.base
        coords = res.coords()
        old_frame = ["C1'"] + list(B.GLYCOSIDIC[old])
        missing = [a for a in old_frame if a not in coords]
        if missing:
            raise B.BuildError(f"{res.label}: missing atom(s) {', '.join(missing)} of the glycosidic frame")
        before, _ = pt.sugar_chirality(coords, old)
        t = B.template(base, hydrogens=True)
        new_frame = ["C1'"] + list(B.GLYCOSIDIC[base])
        mob = np.array([t[a] for a in new_frame])
        ref = np.array([coords[a] for a in old_frame])
        weights = np.array([2.0, 2.0, 1.0, 1.0])
        rot, mc, rc = kabsch(mob, ref, weights)
        rms = float(np.sqrt((((mob - mc) @ rot.T + rc - ref) ** 2).sum(axis=1).mean()))
        template_c1 = next((a for a in res.atoms if a.name.replace("*", "'") == "C1'"), res.atoms[0])
        altlocs = sorted({a.alt for a in res.atoms if a.alt})
        keep, dropped_h, unknown = [], [], []
        insert_at = None
        for a in res.atoms:
            name = a.name.replace("*", "'")
            is_sugar = "'" in name or name in B._SUGAR_EXTRA
            if is_sugar:
                keep.append(a)
            elif name in B._ALL_BASE_NAMES or name in B.BASE_NAMES[old]:
                if insert_at is None:
                    insert_at = len(keep)
                if name.lstrip("0123456789").startswith("H"):
                    dropped_h.append(name)
            else:
                unknown.append(name)
                keep.append(a)
        if insert_at is None:
            insert_at = len(keep)
        record = template_c1.record
        new_names = list(B.BASE_NAMES[base])
        added_h = [h for h, *_ in B._H_RULES[base]] if dropped_h else []
        new_atoms = []
        for name in new_names + added_h:
            xyz = (t[name] - mc) @ rot.T + rc
            tail = template_c1.tail[:3] + B._tail(xyz)[3:33] + template_c1.tail[33:]
            new_atoms.append(B.Atom(record, name, "", tail, B._element(name)))
        res.atoms = keep[:insert_at] + new_atoms + keep[insert_at:]
        old_name = res.name
        res.name = B.mutated_name(res.name, base)
        after, vols = pt.sugar_chirality(res.coords(), base)
        if before is not None and after != before:
            raise B.BuildError(f"{res.label}: the base swap changed the sugar chirality call ({before} -> {after}, "
                               f"{vols})")
        where = f"{res.chain or '_'}:{res.resseq}{res.icode}"
        notes = [f"{where} {old}->{base} ({old_name} -> {res.name}, glycosidic-frame fit {rms:.2f} Å)"]
        if dropped_h:
            notes.append(f"{where}: replaced {len(dropped_h)} hydrogen(s) of the old base ({' '.join(dropped_h)}) "
                         f"with {len(added_h)} of the new base ({' '.join(added_h)}) at ideal positions")
        if unknown:
            notes.append(f"WARNING: {res.label}: kept unrecognised atom(s) {' '.join(unknown)}")
        if altlocs:
            notes.append(f"{where}: alternate locations {','.join(altlocs)}; the new "
                         f"base is fitted on {template_c1.alt or altlocs[0]} and written as one conformer with "
                         "occupancy 1.00")
        return notes, rms

    B._swap_base = _swap_base

    build_helix = B.build_helix

    def build_helix_improved(sequence, partner="auto", chirality="D", chains=("A", "B"), start=(1, 1),
                             scheme=None, terminals=True, hydrogens=False, phosphate5=False, center=True,
                             notes=None, pairing="wc", triplex=None, third="auto"):
        mine = []
        text = build_helix(sequence, partner, chirality, chains, start, scheme, terminals, hydrogens, phosphate5,
                           center, mine, pairing, triplex, third)
        motif = B.check_triplex(triplex)
        seq, other, _ = B.resolve_partner(sequence, partner, B.check_pairing(pairing))
        firsts = B._parse_ids(start, 3 if motif else 2, int)
        used = firsts[:1] + (firsts[1:2] if other is not None else ()) + (firsts[2:3] if motif else ())
        if max(used) + len(seq) - 1 > 9999:  # CH-02
            mine.append("WARNING: residue numbers above 9999 wrap in the PDB columns")
        if notes is not None:
            notes.extend(mine)
        return text.replace("REMARK 951 MDATK IDEAL", "REMARK 951 MDBridge IDEAL")  # D13

    B.build_helix = build_helix_improved

    fit_to_map = B.fit_to_map

    def fit_to_map_checked(source, map_path, resolution, *args, **kwargs):  # CH-10
        if isinstance(resolution, float) and not math.isfinite(resolution):
            raise B.BuildError("the resolution must be a positive number of Å")
        return fit_to_map(source, map_path, resolution, *args, **kwargs)

    B.fit_to_map = fit_to_map_checked


# ---------------------------------------------------------------------------------------------
# the backbone torsion terms of the other pairings' fits (CH-17)

# MDATK's term of a link torsion, (1 + cos 3 theta) * 1.5 + wrap(theta - ref) / 120, jumps by 3 where
# theta - ref crosses +-180 degrees (wrap: +1.5 to -1.5), so the cost is discontinuous there: MDATK's
# reference helices sit next to the jump (hoogsteen strand 2: epsilon +26 / alpha 139; the triplexes'
# alpha 113 and 121, the jump at 118) and closures creep against it until their evaluation limit
# (hoogsteen AAGAGAAGGA: every step at 1/320 of the interval still at its limit; a 2-degree C1 ramp
# in place of the jump left a stiff trough in which hoogsteen GAGGAUCA crept as long). CH-17's tilt is
# smooth and periodic, with MDATK's slope at the reference (1/120 per degree) and MDATK's values at
# the other staggered positions (+-1 at +-120 degrees): TILT_A sin d + TILT_B sin 2d.
TILT_B = (1.5 / math.pi - 2.0 / math.sqrt(3.0)) / 3.0
TILT_A = 2.0 / math.sqrt(3.0) + TILT_B


def torsion_tilt(delta):
    """CH-17's tilt of torsions delta degrees from their references (numpy arrays or floats)."""
    import numpy as np

    r = np.radians(delta)
    return TILT_A * np.sin(r) + TILT_B * np.sin(2.0 * r)


# The reference helices (_helix_fit: every kind's screw and strand torsions) are the one input the
# builds take from them: the uniform helices are them and every closure starts on them. They are
# separable from the closure: _HELIX_START, _SCREW_PRIOR and the torsion term of their fit
# (REFERENCE_TERM: "continuous", CH-17's, or "mdatk", which gives MDATK's reference helices) define
# them; the closure (_patch_closure) always scores the links with CH-17's term. Which reference
# helices are the right ones is open (docs/DIFFERENCES.md CH-17).
REFERENCE_TERM = "continuous"
_MDATK_LINK_RESIDUALS = {}


def continuous_link_residuals(a, b):
    """MDATK's builder._link_residuals with CH-17's tilt."""
    import numpy as np
    from mdatk import builder as B

    def torsion(theta, ref):
        return (1 + np.cos(np.radians(3 * theta))) * 1.5 + torsion_tilt(theta - ref)

    eps = B._dihedral(a["C4'"], a["C3'"], a["O3'"], b["P"])
    return np.concatenate([
        (np.linalg.norm(a["O3'"] - b["P"], axis=-1) - B.LINK) / 0.01,
        (B._angle(a["C3'"], a["O3'"], b["P"]) - B._GEOM["C3'-O3'-P"]) / 1.5,
        (B._angle(a["O3'"], b["P"], b["O5'"]) - B._GEOM["O3'-P-O5'"]) / 1.5,
        torsion(eps, B._TORSIONS["epsilon"]) + np.clip(np.abs(B._wrap(eps + 120.0)) - 60.0, 0.0, None) / 10.0,
        torsion(B._dihedral(a["C3'"], a["O3'"], b["P"], b["O5'"]), B._TORSIONS["zeta"]),
        torsion(B._dihedral(a["O3'"], b["P"], b["O5'"], b["C5'"]), B._TORSIONS["alpha"])])


def _patch_torsions():
    """The reference fit's term (builder._link_residuals, which _fit_residuals calls) as REFERENCE_TERM."""
    from mdatk import builder as B

    _MDATK_LINK_RESIDUALS.setdefault("f", B._link_residuals)
    B._link_residuals = continuous_link_residuals if REFERENCE_TERM == "continuous" else _MDATK_LINK_RESIDUALS["f"]


# ---------------------------------------------------------------------------------------------
# the other pairings' strand closure (CH-16, owner decision D17)

# Continuation: the partner nucleosides move from their strand's reference pose (t = 0, the uniform
# helix closed by _helix_fit's torsions) to their _level_pose poses (t = 1) along the rigid transform
# (rotation vector about the reference C1' and C1' shift, both scaled by t). The interval is walked in
# CLOSURE_STEPS equal steps; a step whose solve stops on its evaluation limit, or whose contact pairs
# needed more than CLOSURE_ROUNDS_OK refreshes, is retried as two half steps (from the last accepted
# point, recursively, down to CLOSURE_MAX_HALVINGS halvings; the smallest steps get the larger limit
# and are accepted as they end). Every solve: MDATK's residuals at its last pass's contact weight, TRF
# with the exact subproblem (2-point differences grouped by MDATK's level-span sparsity, then dense),
# ftol = xtol = gtol = CLOSURE_TOL; the contact pairs within CLOSURE_CUTOFF (MDATK's filters) are
# refreshed before every solve, kept once scored, in (p, q) order, and the solve repeats while the set
# grows (at most CLOSURE_MAX_ROUNDS times). A solve's variables are the move from its start (TRF's
# first trust radius is then 1 in x_scale's units, not the norm of the start: a start of torsions near
# 100 degrees let the first steps jump by hundreds of degrees into another basin).
CLOSURE_STEPS = 10
CLOSURE_MAX_HALVINGS = 5          # down to 1 / 320 of the interval
CLOSURE_WEIGHT = 0.01             # contact weight (MDATK's last pass)
CLOSURE_TOL = 1e-10
CLOSURE_MAX_NFEV = 1000           # evaluation limit of a step that may still be halved
CLOSURE_FINAL_NFEV = 20000        # ... of the t = 0 solve and of the smallest steps
CLOSURE_CUTOFF = 4.5              # Å (MDATK: 3.8)
CLOSURE_ROUNDS_OK = 4             # pair refreshes a step may need before it is halved
CLOSURE_MAX_ROUNDS = 50
SMALL_FIT_TOL = 1e-14             # _pair_fit, _level_pose, _helix_fit
SMALL_FIT_MAX_NFEV = 5000

# the solver the closure and the small fits call (a test harness may wrap it)
SOLVER = {}


def _rotation_vector(R):
    """Rotation vector (radians) of a rotation matrix (angle below pi)."""
    import numpy as np

    c = (np.trace(R) - 1.0) / 2.0
    th = math.acos(min(1.0, max(-1.0, c)))
    if th < 1e-12:
        return np.zeros(3)
    w = np.array([R[2, 1] - R[1, 2], R[0, 2] - R[2, 0], R[1, 0] - R[0, 1]]) / (2.0 * math.sin(th))
    return w * th


def column_groups(structure):
    """Columns of a boolean Jacobian structure (m, n) that share no row, as SciPy's
    group_columns(structure, order=0) forms them (greedy over numpy.random.RandomState(0).permutation(n))."""
    import numpy as np

    m, n = structure.shape
    order = np.random.RandomState(0).permutation(n)
    cols = structure[:, order].T  # (n, m)
    groups = -np.ones(n, dtype=np.intp)
    current = 0
    for i in range(n):
        if groups[i] >= 0:
            continue
        groups[i] = current
        union = cols[i].copy()
        all_grouped = True
        for j in range(n):
            if groups[j] >= 0:
                continue
            all_grouped = False
            if not np.any(union & cols[j]):
                union |= cols[j]
                groups[j] = current
        if all_grouped:
            break
        current += 1
    out = np.empty(n, dtype=np.intp)
    out[order] = groups
    return out


def fd_jacobian(fun, x, rows_of, groups):
    """2-point finite-difference Jacobian (dense) with column groups (rows_of[j]: the rows of column j in
    the structure), as SciPy's approx_derivative(fun, x, method="2-point", sparsity=(structure, groups))
    computes it."""
    import numpy as np

    f0 = fun(x)
    h = np.finfo(np.float64).eps ** 0.5 * np.where(x >= 0, 1.0, -1.0) * np.maximum(1.0, np.abs(x))
    jac = np.zeros((len(f0), len(x)))
    for g in range(int(groups.max()) + 1):
        cols = np.nonzero(groups == g)[0]
        xg = x.copy()
        xg[cols] = x[cols] + h[cols]
        df = fun(xg) - f0
        for j in cols:
            rows = rows_of[j]
            jac[rows, j] = df[rows] / (xg[j] - x[j])
    return jac


def _patch_closure():
    import functools

    import numpy as np
    import scipy.optimize
    from mdatk import builder as B
    from mdatk.spatial import kabsch

    SOLVER.setdefault("least_squares", scipy.optimize.least_squares)
    depth = [0]

    def least_squares(fun, x0, *args, **kwargs):  # the small fits converged
        if depth[0] > 0:
            kwargs.update(ftol=SMALL_FIT_TOL, xtol=SMALL_FIT_TOL, gtol=SMALL_FIT_TOL,
                          max_nfev=kwargs.get("max_nfev") or SMALL_FIT_MAX_NFEV)
        return SOLVER["least_squares"](fun, x0, *args, **kwargs)

    scipy.optimize.least_squares = least_squares  # the builder imports it when it fits

    def converged(fn):
        inner = getattr(fn, "__wrapped__", fn)

        @functools.lru_cache(maxsize=None)
        def wrapper(*args):
            depth[0] += 1
            try:
                return inner(*args)
            finally:
                depth[0] -= 1
        return wrapper

    for name in ("_pair_fit", "_level_pose", "_helix_fit"):
        setattr(B, name, converged(getattr(B, name)))

    context = {}
    helix_strands = B._helix_strands

    def _helix_strands(kind, sequences):
        context["kind"] = kind
        return helix_strands(kind, sequences)

    B._helix_strands = _helix_strands
    names, slots, base_code = B._STRAND_NAMES, B._BASE_SLOTS, B._BASE_CODE

    def _close_levels(placed, bases, steps, start):
        kind = context["kind"]
        spec = B._HELICES[kind]
        screw = B._helix_fit(kind).screw
        n_s, n = len(placed), len(placed[0])
        width = 6 + 3 * n_s
        # per partner nucleoside: the reference pose and the rigid move to its _level_pose pose
        moves = [None]
        for s in range(1, n_s):
            rot, shift = B._pair_fit(spec[s].rule, spec[0].base, spec[s].base)
            row = []
            for k in range(n):
                ref = {a: screw.apply(v, k) for a, v in B._moved(B._core(bases[s][k]), rot, shift).items()}
                common = [a for a in placed[s][k] if a in ref and a != "H3+"]
                P = np.array([ref[a] for a in common])
                Q = np.array([placed[s][k][a] for a in common])
                R, pc, qc = kabsch(P, Q, np.ones(len(common)))
                c = ref["C1'"]
                row.append((ref, _rotation_vector(R), c, (c - pc) @ R.T + qc - c))
            moves.append(row)

        def placed_at(t):
            if t == 1.0:
                return placed
            out = [placed[0]]
            for s in range(1, n_s):
                row = []
                for ref, w, c, d in moves[s]:
                    R = B._rotvec(w * t)
                    row.append({a: (v - c) @ R.T + c + t * d for a, v in ref.items()})
                out.append(row)
            return out

        codes = np.array([[base_code[b] for b in bs] for bs in bases])
        cls_tab, inside_tab, link_tab = B._contact_tables()
        n_atoms = len(names) + slots
        cls = cls_tab[codes].ravel()
        atom_s, atom_k, atom_i = np.unravel_index(np.arange(n_s * n * n_atoms), (n_s, n, n_atoms))
        step_of = np.asarray(steps)
        order = [np.arange(n) if st > 0 else np.arange(n)[::-1] for st in steps]
        geo = {}

        def set_geometry(pl):
            arrays = [B._core_arrays(p, b, True) for p, b in zip(pl, bases)]
            base_xyz = np.full((n_s, n, slots, 3), 1e4)
            for s in range(n_s):
                for k, (core, b) in enumerate(zip(pl[s], bases[s])):
                    base_xyz[s, k, :len(B._base_slots(b))] = [core[a] for a in B._base_slots(b)]
            geo.update(arrays=arrays, base_xyz=base_xyz, centre=np.mean([a[0] for a in arrays], axis=0))

        def build(x):  # MDATK's
            v = x.reshape(n, width)
            rot = B._rotvecs(np.radians(v[:, 0:3]))
            centre = geo["centre"]

            def move(y):
                return np.einsum("nij,nmj->nmi", rot, y - centre[:, None]) + centre[:, None] + v[:, None, 3:6]

            nts, atoms = [], []
            for s, (c1, glyco, sugar) in enumerate(geo["arrays"]):
                q = v[:, 6 + 3 * s:9 + 3 * s]
                nt = B._nucleotides(move(c1[:, None])[:, 0], move(glyco[:, None])[:, 0], move(sugar), q[:, 0],
                                    q[:, 1], q[:, 2])
                o = order[s]
                prev = np.empty((n, 3))
                prev[o[1:]] = nt["O3'"][o[:-1]]
                prev[o[0]] = B._place(nt["C5'"][o[0]], nt["O5'"][o[0]], nt["P"][o[0]], B.LINK,
                                      B._GEOM["O3'-P-O5'"], B._TORSIONS["alpha"])
                nt["OP1"], nt["OP2"] = B._phosphate_oxygens_many(nt["P"], prev, nt["O5'"])
                nts.append(nt)
                atoms.append(np.concatenate([np.stack([nt[a] for a in names], axis=1), move(geo["base_xyz"][s])],
                                            axis=1))
            return nts, np.array(atoms).reshape(-1, 3)

        pairs = {}

        def refresh(x):
            """The pairs within CLOSURE_CUTOFF (MDATK's filters) added to those scored so far, in (p, q)
            order; whether the set grew."""
            from scipy.spatial import cKDTree

            xyz = build(x)[1]
            real = np.nonzero(cls != 6)[0]
            found = cKDTree(xyz[real]).query_pairs(CLOSURE_CUTOFF, output_type="ndarray")
            p, q = real[found[:, 0]], real[found[:, 1]]
            sp, kp, ip, sq, kq, iq = atom_s[p], atom_k[p], atom_i[p], atom_s[q], atom_k[q], atom_i[q]
            keep = (kp != kq) | (sp != sq) & ((ip < len(names)) | (iq < len(names)))
            same = (sp == sq) & (kp == kq)
            keep[same] = inside_tab[codes[sp[same], kp[same]], ip[same], iq[same]]
            linked = (sp == sq) & (kq - kp == step_of[sp])
            keep[linked] = link_tab[codes[sp[linked], kp[linked]], codes[sq[linked], kq[linked]], ip[linked],
                                    iq[linked]]
            linked = (sp == sq) & (kp - kq == step_of[sp])
            keep[linked] = link_tab[codes[sq[linked], kq[linked]], codes[sp[linked], kp[linked]], iq[linked],
                                    ip[linked]]
            limit = B._FIT_LIMITS[cls[p], cls[q]]
            keep &= limit > 0
            ext_p = np.clip(np.where((ip == 1) | (ip == 2), kp - step_of[sp], kp), 0, n - 1)
            ext_q = np.clip(np.where((iq == 1) | (iq == 2), kq - step_of[sq], kq), 0, n - 1)
            lo = np.minimum(np.minimum(kp, kq), np.minimum(ext_p, ext_q))
            hi = np.maximum(np.maximum(kp, kq), np.maximum(ext_p, ext_q))
            new = dict(p=p[keep], q=q[keep], limit=limit[keep], lo=lo[keep], hi=hi[keep])
            before = len(pairs.get("p", ()))
            if pairs:
                new = {k: np.concatenate([pairs[k], new[k]]) for k in new}
            _u, first = np.unique(new["p"] * len(cls) + new["q"], return_index=True)
            pairs.update({k: v[first] for k, v in new.items()})
            return len(pairs["p"]) != before

        def residuals(x):  # MDATK's
            v = x.reshape(n, width)
            nts, xyz = build(x)
            out = [continuous_link_residuals({k: a[o[:-1]] for k, a in nt.items()},
                                             {k: a[o[1:]] for k, a in nt.items()})
                   for nt, o in zip(nts, order)]
            out += [B._torsion_priors(v[:, 6 + 3 * s], v[:, 7 + 3 * s], v[:, 8 + 3 * s]) for s in range(n_s)]
            viol = pairs["limit"] - np.linalg.norm(xyz[pairs["p"]] - xyz[pairs["q"]], axis=-1)
            out += [(v[:, 0:3] / B.LEVEL_MOVE[0]).ravel(), (v[:, 3:6] / B.LEVEL_MOVE[1]).ravel(),
                    np.clip(viol, 0.0, None) / CLOSURE_WEIGHT]
            return np.concatenate(out)

        def structure():  # MDATK's level spans
            lo, hi = [], []
            for o in order:
                a, b = np.minimum(o[:-1], o[1:]), np.maximum(o[:-1], o[1:])
                lo += list(a) * 6
                hi += list(b) * 6
            lo += list(range(n)) * 4 * n_s
            hi += list(range(n)) * 4 * n_s
            lo += list(np.repeat(np.arange(n), 3)) * 2
            hi += list(np.repeat(np.arange(n), 3)) * 2
            lo = np.concatenate([np.asarray(lo, dtype=np.intp), pairs["lo"]])
            hi = np.concatenate([np.asarray(hi, dtype=np.intp), pairs["hi"]])
            level = np.arange(n * width) // width
            return (lo[:, None] <= level[None]) & (level[None] <= hi[:, None])

        scale = np.tile([3.0, 3.0, 3.0, 0.2, 0.2, 0.2] + [10.0] * 3 * n_s, n)
        group_cache = {}

        def solve(x, max_nfev):
            sp = structure()
            key = sp.tobytes()
            if key not in group_cache:
                group_cache.clear()
                group_cache[key] = (column_groups(sp), [np.nonzero(sp[:, j])[0] for j in range(sp.shape[1])])
            groups, rows_of = group_cache[key]

            def moved(y):  # the variables as the move from x: TRF's first trust radius is then 1 (scaled)
                return residuals(x + y)

            sol = SOLVER["least_squares"](moved, np.zeros_like(x),
                                          jac=lambda yy: fd_jacobian(moved, yy, rows_of, groups),
                                          tr_solver="exact", x_scale=scale, ftol=CLOSURE_TOL, xtol=CLOSURE_TOL,
                                          gtol=CLOSURE_TOL, max_nfev=max_nfev)
            return x + sol.x, sol.status == 0

        def step(x, t, max_nfev, halvable):
            """Solve at t from x with the pair refreshes: (x, whether the step had trouble); a step that
            can still be halved stops at its first trouble."""
            set_geometry(placed_at(t))
            trouble = False
            for rounds in range(CLOSURE_MAX_ROUNDS):
                if not refresh(x) and rounds > 0:
                    break
                x, hit = solve(x, max_nfev)
                trouble = trouble or hit or rounds >= CLOSURE_ROUNDS_OK
                if trouble and halvable:
                    break
            return x, trouble

        def advance(x, t0, t1, halvings):
            last = halvings >= CLOSURE_MAX_HALVINGS
            saved = {k: v.copy() for k, v in pairs.items()}
            y, trouble = step(x, t1, CLOSURE_FINAL_NFEV if last else CLOSURE_MAX_NFEV, not last)
            if not trouble or last:
                return y
            pairs.clear()
            pairs.update(saved)
            mid = 0.5 * (t0 + t1)
            return advance(advance(x, t0, mid, halvings + 1), mid, t1, halvings + 1)

        x = np.concatenate([np.zeros((n, 6))] + [np.asarray(q, dtype=np.float64).reshape(n, 3) for q in start],
                           axis=1).ravel()
        x = step(x, 0.0, CLOSURE_FINAL_NFEV, False)[0]
        for i in range(CLOSURE_STEPS):
            x = advance(x, i / CLOSURE_STEPS, (i + 1) / CLOSURE_STEPS, 0)
        set_geometry(placed)
        v = x.reshape(n, width)
        rot = B._rotvecs(np.radians(v[:, 0:3]))
        centre = geo["centre"]
        moved = [[{a: (y - centre[k]) @ rot[k].T + centre[k] + v[k, 3:6] for a, y in core.items()}
                  for k, core in enumerate(strand)] for strand in placed]
        return [v[:, 6 + 3 * s:9 + 3 * s] for s in range(n_s)], moved

    B._close_levels = _close_levels


# ---------------------------------------------------------------------------------------------
# density (CH-11)

def finite_map(dmap):
    """(the map with its non-finite voxels as 0, their count)."""
    import numpy as np
    from mdatk.io.mrc import DensityMap

    bad = ~np.isfinite(dmap.data)
    n = int(bad.sum())
    if not n:
        return dmap, 0
    return DensityMap(np.where(bad, 0, dmap.data).astype(dmap.data.dtype), dmap.origin, dmap.voxel), n


def finite_level(dmap, sigmas):
    """mean + sigmas * rms of the finite voxels (RN-33)."""
    import numpy as np

    v = dmap.data[np.isfinite(dmap.data)].astype(np.float64)
    if not len(v):
        return float("nan")
    return float(v.mean() + sigmas * v.std())


def _patch_density():
    from mdatk import density as D

    fit_model = D.fit_model

    def fit_model_finite(dmap, xyz, elements, resolution, start=None, search="global", angular_step=30.0,
                         progress=None, cancelled=None):
        dmap, n = finite_map(dmap)
        res = fit_model(dmap, xyz, elements, resolution, start=start, search=search, angular_step=angular_step,
                        progress=progress, cancelled=cancelled)
        if n:
            res.message = (res.message + "; " if res.message else "") + \
                f"{n} non-finite map voxel(s) read as no density"
        return res

    D.fit_model = fit_model_finite
    D._mdbridge_finite = True


# ---------------------------------------------------------------------------------------------
# the command line

def _missing_input(path):
    """CH-13: OSError for a path that is not a file."""
    if path is None or os.path.isfile(path):
        return
    if os.path.isdir(path):
        raise IsADirectoryError(errno.EISDIR, os.strerror(errno.EISDIR), path)
    raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), path)


def apply_cli():
    global _applied_cli
    if _applied_cli:
        return
    _applied_cli = True
    apply_mdatk()
    import numpy as np  # noqa: F401
    from mdatk import builder as B
    from mdatk import cli, dictionaries, pipeline
    from mdatk import pdbtools as pt

    pdb_run = pt.run

    def pdb_run_checked(args, stdout=None, stderr=None):  # CH-13
        _missing_input(args.input)
        if getattr(args, "chain_reference", None):
            _missing_input(args.chain_reference)
        return pdb_run(args, stdout, stderr)

    pt.run = pdb_run_checked

    annotate = cli._pdb_annotate

    def pdb_annotate_checked(args, name):
        _missing_input(args.input)
        return annotate(args, name)

    cli._pdb_annotate = pdb_annotate_checked

    build_run = B.run

    def build_run_checked(args, progress=None):
        if args.build_command in ("mutate", "fit"):
            _missing_input(args.input)
        return build_run(args, progress)

    B.run = build_run_checked

    # D7: the PubChem root of the stand-in server
    fetch = dictionaries.fetch_pubchem_smiles

    def fetch_pubchem_smiles(cid_or_name, timeout=15.0):
        url = os.environ.get(PUBCHEM_URL_VARIABLE)
        saved = dictionaries.PUBCHEM_REST
        if url:
            dictionaries.PUBCHEM_REST = url.rstrip("/")
        try:
            return fetch(cid_or_name, timeout)
        finally:
            dictionaries.PUBCHEM_REST = saved

    dictionaries.fetch_pubchem_smiles = fetch_pubchem_smiles

    # CH-11: the density centre of a local start from the finite voxels
    fit_density = pipeline.Pipeline.fit_density

    def fit_density_finite(self, search=None, start=None, reporter=None):
        from mdatk import density as D
        level, centre = D.map_level, D.centre_transform

        def centre_transform(dmap, lvl, model_centre):
            return centre(finite_map(dmap)[0], lvl, model_centre)

        D.map_level, D.centre_transform = finite_level, centre_transform
        try:
            return fit_density(self, search=search, start=start, reporter=reporter)
        finally:
            D.map_level, D.centre_transform = level, centre

    pipeline.Pipeline.fit_density = fit_density_finite

    cli.cmd_density_fit = _cmd_density_fit
    simulate = cli.cmd_density_simulate

    def cmd_density_simulate(args):  # CH-10
        if not (math.isfinite(args.resolution) and args.resolution > 0):
            raise ValueError("the resolution must be a positive number of Å")
        if args.spacing is not None and not (math.isfinite(args.spacing) and args.spacing >= 0):
            raise ValueError("the grid spacing must be positive")
        return simulate(args)

    cli.cmd_density_simulate = cmd_density_simulate  # (build_parser reads the module globals)


def _cmd_density_fit(args):
    """MDATK 0.0.4.1's cmd_density_fit with CH-10 and CH-12."""
    import json

    from mdatk import cli
    resolution = getattr(args, "resolution", None)
    if resolution is not None and not math.isfinite(resolution):  # CH-10
        raise ValueError("the resolution must be a positive number of Å")
    project, pipe = cli._pipeline(args)
    cli._density_overrides(project, args)
    problems = project.density_problems()
    if not project.density.map_path:
        problems.append("no density map: give --map or set density.map_path in the project")
    if problems:
        raise ValueError("; ".join(problems))
    pipe.check()
    res = pipe.fit_density(search=args.search)
    d = project.density
    out = {"cc": round(float(res.cc), 4), "start": res.start, "info": d.fit_info, "transform": d.transform}
    if args.save:
        target = args.save if isinstance(args.save, str) else args.params
        if not target:
            raise ValueError("--save needs a path when no --project was given")
        if args.params:
            _save_as_written(args, project, target)
        else:
            project.save(target)
        out["saved"] = os.path.abspath(target)
    print(json.dumps(out, indent=2))
    return 0


def _save_as_written(args, project, target):
    """CH-12: the project file as written, with the fit's density settings (ED-02; Save As rebases)."""
    import improvements_gui
    from mdatk import projectfile
    from mdatk.config import Project

    source = os.path.abspath(args.params)
    written = improvements_gui.load_as_written(Project, projectfile, source)
    d, w = project.density, written.density
    if args.map:
        given = args.map
        w.map_path = given if os.path.isabs(given) else os.path.relpath(os.path.abspath(given),
                                                                         os.path.dirname(source))
    for name in ("resolution", "level", "fit_frame", "angular_step", "transform", "anchor", "fit_cc", "fit_info"):
        setattr(w, name, copy.deepcopy(getattr(d, name)))
    old_dir, new_dir = os.path.dirname(source), os.path.dirname(os.path.abspath(target))
    if old_dir != new_dir:
        written = improvements_gui.rebase(written, old_dir, new_dir)
    written.save(target)
