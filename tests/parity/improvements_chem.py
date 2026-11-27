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
  CH-17  the link torsions' term of the other pairings' closure is continuous (continuous_link_residuals:
         MDATK's (1 + cos 3 theta) * 1.5 + wrap(theta - ref) / 120 jumps by 3 at theta = ref + 180,
         where closures crept): the tilt TILT_A sin d + TILT_B sin 2d (d = theta - ref) has MDATK's
         slope at the reference and its values +-1 at the staggered positions +-120 degrees; the
         reference helices' fits keep MDATK's term (D18); 0.1.1 has MDATK's term
  CH-18  (owner decision D18) realistic non-Watson-Crick geometry (_patch_ch18): RNA builds of
         reverse-wc, reverse-hoogsteen and the reverse-hoogsteen (purine-motif) triplex are refused
         (refusal(): build_helix raises, check_helix gives the one error, after MDATK's own request
         checks); the Hoogsteen duplex and the pyrimidine-motif triplex take their nucleotides, pairs,
         screw and stacking from the RNA triple helices of PDB 3P22, 4PLX, 6SVS, 7JNH
         (tools/gen_chem_tables.triplex_tables; strand nucleotides exp_core, pairs table_pose, the
         restrained reference fit exp_fit), and their closures keep the stacking; every reference fit
         scores contacts over CONTACT_LEVELS levels (MDATK: 2); a build with a partner notes the pairs
         without a base-base H-bond (hbond_note); REMARK 951 / the build note name the triplex
         structures; 0.1.1 has MDATK's builder
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
    _patch_closure()
    _patch_ch18()


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

    helix_strands = B._helix_strands

    def _helix_strands(kind, sequences):
        CLOSURE_CONTEXT["kind"] = kind
        return helix_strands(kind, sequences)

    B._helix_strands = _helix_strands
    names, slots, base_code = B._STRAND_NAMES, B._BASE_SLOTS, B._BASE_CODE

    def _close_levels(placed, bases, steps, start):
        kind = CLOSURE_CONTEXT["kind"]
        spec = B._HELICES[kind]
        screws = triplex()["screws"]
        stack = (screws[kind]["stack"], screws[kind]["stack_sd"]) if kind in EXP_KINDS else None
        screw = B._helix_fit(kind).screw
        n_s, n = len(placed), len(placed[0])
        width = 6 + 3 * n_s
        # per partner nucleoside: the reference pose and the rigid move to its _level_pose pose
        moves = [None]
        for s in range(1, n_s):
            rot, shift = reference_pose(kind, s)  # CH-18: the kind's reference pose and nucleosides
            row = []
            for k in range(n):
                ref = {a: screw.apply(v, k) for a, v in B._moved(kind_core(kind, s, bases[s][k]), rot, shift).items()}
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
            out += [(v[:, 0:3] / B.LEVEL_MOVE[0]).ravel(), (v[:, 3:6] / B.LEVEL_MOVE[1]).ravel()]
            if stack:  # CH-18: the base-plane separation of consecutive levels
                a = xyz.reshape(n_s, n, n_atoms, 3)
                levels = [np.concatenate([a[s_, k, len(names):len(names) + len(B.BASE_NAMES[bases[s_][k]])]
                                          for s_ in range(n_s)]) for k in range(n)]
                out.append(np.array([(plane_separation(levels[k], levels[k + 1]) - stack[0]) / stack[1]
                                     for k in range(n - 1)]))
            out.append(np.clip(viol, 0.0, None) / CLOSURE_WEIGHT)
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
            if stack:
                lo += list(range(n - 1))
                hi += list(range(1, n))
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
# realistic non-Watson-Crick geometry (CH-18, owner decision D18)

# The kinds RNA is not known to form are refused: what the error names, per pairing / triplex motif.
REFUSED_PAIRINGS = {"reverse-wc": "the parallel reverse Watson-Crick (tWW) duplex",
                    "reverse-hoogsteen": "the reverse-Hoogsteen (tHW) duplex"}
REFUSED_TRIPLEXES = {"reverse-hoogsteen": "the purine-motif (reverse-Hoogsteen) triplex"}
# The kinds whose reference geometry comes from the RNA triple helices (gen_chem_tables.triplex_tables):
# kind -> the table strand (0 purine, 1 Watson-Crick, 2 Hoogsteen) of each of its strands.
EXP_KINDS = {"hoogsteen": (0, 2), "triplex-hoogsteen": (0, 1, 2)}
# the canonical pairs the tables give: (rule, strand-1 base, partner) -> (level type, table strand)
TABLE_PAIRS = {("wc", "A", "U"): ("AUU", 1), ("hoogsteen", "A", "U"): ("AUU", 2),
               ("wc", "G", "C"): ("GCC", 1), ("hoogsteen", "G", "C"): ("GCC", 2)}
GC_HBOND = (2.8, 3.1, 0.05)  # Å: the band the refined C+·G-C pairs' H-bonds are brought into, the scale of a miss
CONTACT_LEVELS = 12  # levels after a reference level whose contacts the reference fits score (a turn)
REFERENCE_CUTOFF = 6.0  # Å: the contact pairs a restrained reference fit scores (refreshed at its solutions)
REFERENCE_ROUNDS = 20
REFERENCE_STARTS = ((0.0, 58.0, 174.0), (-10.0, 60.0, 160.0), (10.0, 50.0, 190.0), (0.0, 80.0, 150.0),
                    (0.0, 120.0, 180.0))  # torsion starts (chi change, gamma, beta) of the restrained fits
REFERENCE_SCALE_SCREW = (1.0, 1.0, 5.0, 5.0, 5.0, 0.3)
REFERENCE_SCALE_POSE = (3.0, 3.0, 3.0, 0.2, 0.2, 0.2)
HBOND_NOTE = 3.5  # Å: a base-base H-bond of the note is a donor-acceptor N/O pair closer than this
EXP_SOURCES = "3P22, 4PLX, 6SVS and 7JNH"


def refusal(pairing, motif):
    """The error of a kind RNA is not known to form, or None."""
    what = REFUSED_TRIPLEXES.get(motif) if motif else REFUSED_PAIRINGS.get(pairing)
    if what is None:
        return None
    asked = f"triplex {motif}" if motif else f"pairing {pairing}"
    return (f"{asked}: RNA is not known to form {what}; DNA and DNA-RNA hybrid builds of it will come with "
            "DNA structure building")


def triplex():
    if "t" not in _TRIPLEX:
        _TRIPLEX["t"] = chem_tables().triplex_tables()
    return _TRIPLEX["t"]


_TRIPLEX = {}
_CACHE = {}


def _cached(key, make):
    if key not in _CACHE:
        _CACHE[key] = make()
    return _CACHE[key]


def _table_atoms(level_type, strand):
    import numpy as np

    return {n: np.array([x, y, z]) for n, x, y, z in triplex()["levels"][level_type][strand]}


def strand_nucleotide(t):
    """The template nucleotide of table strand t in the helix frame of the tables: the 4JRD nucleotide
    (CHD-02) without its base, superposed by C1' and its sugar atoms onto the mean U·A-U nucleoside of
    the strand, and that nucleoside's mean base atoms."""
    import numpy as np
    from mdatk import builder as B
    from mdatk.spatial import kabsch

    def make():
        base = "AUU"[t]
        table = _table_atoms("AUU", t)
        jrd = {n: B._cartesian(r, p, z) for n, r, p, z in B._JRD_NUCLEOTIDE}
        names = ["C1'"] + list(chem_tables().TRIPLEX_SUGAR)
        rot, mc, rc = kabsch(np.array([jrd[n] for n in names]), np.array([table[n] for n in names]),
                             np.ones(len(names)))
        out = {n: (x - mc) @ rot.T + rc for n, x in jrd.items() if n not in B.BASE_NAMES["A"]}
        out.update({n: table[n] for n in B.BASE_NAMES[base]})
        return base, out
    return _cached(("nucleotide", t), make)


def exp_core(t, base):
    """_core_template on the strand-t template nucleotide (the Arnott base of `base` on its glycosidic
    frame; base and ring hydrogens; C: H3+)."""
    import numpy as np
    from mdatk import builder as B
    from mdatk.spatial import kabsch

    def make():
        tb, nuc = strand_nucleotide(t)
        arn = B.template(base, hydrogens=False)
        old, new = ["C1'"] + list(B.GLYCOSIDIC[tb]), ["C1'"] + list(B.GLYCOSIDIC[base])
        rot, mc, rc = kabsch(np.array([arn[a] for a in new]), np.array([nuc[a] for a in old]),
                             np.array([2.0, 2.0, 1.0, 1.0]))
        heavy = {n: x for n, x in nuc.items() if n not in B.BASE_NAMES[tb]}
        heavy.update({n: (arn[n] - mc) @ rot.T + rc for n in B.BASE_NAMES[base]})
        keep = {h for h, *_ in B._H_RULES[base]} | set(B._SUGAR_H)
        hyd = {n: x for n, x in B._hydrogens(heavy, base).items() if n in keep}
        if base == "C":
            n3 = heavy["N3"]
            hyd["H3+"] = n3 - B.BOND_H["nh"] * B._unit(B._unit(heavy["C2"] - n3) + B._unit(heavy["C4"] - n3))
        return tuple((n, tuple(x)) for n, x in list(heavy.items()) + list(hyd.items()))
    return {n: np.array(x) for n, x in _cached(("core", t, base), make)}


def kind_core(kind, strand, base):
    from mdatk import builder as B

    return exp_core(EXP_KINDS[kind][strand], base) if kind in EXP_KINDS else B._core(base)


def table_pose(rule, b1, b2):
    """(rotation, shift) of the strand core of b2 on the strand-1 core of b1 as in the mean level of
    the tables (each core superposed, unweighted, on C1', its base and sugar atoms there)."""
    import numpy as np
    from mdatk import builder as B
    from mdatk.spatial import kabsch

    def make():
        level, t = TABLE_PAIRS[(rule, b1, b2)]
        names1, names2 = chem_tables().nucleoside_names(b1), chem_tables().nucleoside_names(b2)
        a1, a2 = _table_atoms(level, 0), _table_atoms(level, t)
        c1, c2 = exp_core(0, b1), exp_core(t, b2)
        r1, m1, s1 = kabsch(np.array([c1[n] for n in names1]), np.array([a1[n] for n in names1]),
                            np.ones(len(names1)))
        r2, m2, s2 = kabsch(np.array([c2[n] for n in names2]), np.array([a2[n] for n in names2]),
                            np.ones(len(names2)))
        rot0 = r1.T @ r2
        shift0 = (s2 - r2 @ m2 - s1) @ r1 + m1
        if level == "AUU":
            return rot0, shift0
        # the C+·G-C mean (5 triples; its N2-O2 2.58 Å, its O6-H41 1.6 Å) refined: the partner moved
        # (about its C1') as little as the U·A-U pose spreads allow while its H-bonds come into
        # the band GC_HBOND
        rot_sd, shift_sd = {st: (r, sh) for st, r, sh in triplex()["pose_sd"]}[t]
        start = {n: x @ rot0.T + shift0 for n, x in c2.items()}
        centre = start["C1'"]
        bonds = B._PAIR_HBONDS[(rule, b1, b2)]

        def residuals(p):
            r = B._rotvec(np.radians(p[:3]))
            res = list(p[:3] / rot_sd) + list(p[3:] / shift_sd)
            for x, y, *_ in bonds:
                d = np.linalg.norm(c1[x] - ((start[y] - centre) @ r.T + centre + p[3:]))
                res.append((max(GC_HBOND[0] - d, 0.0) + max(d - GC_HBOND[1], 0.0)) / GC_HBOND[2])
            return np.array(res)

        sol = SOLVER["least_squares"](residuals, np.zeros(6), ftol=SMALL_FIT_TOL, xtol=SMALL_FIT_TOL,
                                      gtol=SMALL_FIT_TOL, max_nfev=SMALL_FIT_MAX_NFEV)
        r = B._rotvec(np.radians(sol.x[:3]))
        return r @ rot0, (shift0 - centre) @ r.T + centre + sol.x[3:]
    return _cached(("table pose", rule, b1, b2), make)


def hbond_pair_fit(rule, b1, b2, t1, t2):
    """MDATK's _pair_fit on the given nucleosides t1 (strand 1) and t2 (its body, MDATK 0.0.4.1)."""
    import numpy as np
    from scipy.optimize import least_squares
    from mdatk import builder as B
    from mdatk.spatial import kabsch

    bonds = B._PAIR_HBONDS[(rule, b1, b2)]
    ring1, n1 = B._base_plane(t1, b1)
    ring2, n2 = B._base_plane(t2, b2)
    a = np.array([t1[x] for x, *_ in bonds])
    b = np.array([t2[y] for _x, y, *_ in bonds])
    out = B._unit(np.cross(n1, a[-1] - a[0]))
    if np.dot(out, a.mean(axis=0) - ring1) < 0:
        out = -out
    targets = a + B.HBOND * out
    centre = targets.mean(axis=0) + out * np.linalg.norm(b.mean(axis=0) - ring2)
    names1 = list(B.BASE_NAMES[b1]) + ["C1'"] + list(B._SUGAR)
    names2 = list(B.BASE_NAMES[b2]) + ["C1'"] + list(B._SUGAR)
    fixed = np.array([t1[n] for n in names1])
    mobile = np.array([t2[n] for n in names2])
    hbonded = {(x, y) for x, y, *_ in bonds}
    free = np.array([[(x, y) not in hbonded for y in names2] for x in names1])
    planar = np.array([t2[n] for n in B.BASE_NAMES[b2]])
    best = None
    for sign in (1.0, -1.0):
        r0, cs, cd = kabsch(np.vstack([b, ring2, ring2 + n2]), np.vstack([targets, centre, centre + sign * n1]),
                            np.array([10.0] * len(b) + [1.0, 1.0]))

        def move(p, x, r0=r0, cs=cs, cd=cd):
            return (x - cs) @ (B._rotvec(p[:3]) @ r0).T + cd + p[3:]

        def residuals(p, move=move):
            res = []
            for x, y, side, h in bonds:
                acc, don = (move(p, t2[y]), t1[x]) if side == 1 else (t1[x], move(p, t2[y]))
                hyd = t1[h] if side == 1 else move(p, t2[h])
                d = np.linalg.norm(don - acc)
                res += [(d - B.HBOND) * 4.0, (d - np.linalg.norm(hyd - acc) - B.BOND_H["nh"]) * 3.0]
            res += list((move(p, planar) - ring1) @ n1)
            dist = np.linalg.norm(fixed[:, None] - move(p, mobile)[None], axis=-1)
            res += list((np.clip(3.1 - dist, 0.0, None) * free).ravel() * 3.0)
            return np.array(res)

        sol = SOLVER["least_squares"](residuals, np.zeros(6), ftol=SMALL_FIT_TOL, xtol=SMALL_FIT_TOL,
                                      gtol=SMALL_FIT_TOL, max_nfev=SMALL_FIT_MAX_NFEV)
        if best is None or sol.cost < best[0]:
            rot = B._rotvec(sol.x[:3]) @ r0
            best = (sol.cost, rot, cd + sol.x[3:] - cs @ rot.T)
    return best[1], best[2]


def kind_pair_fit(kind, strand, rule, b1, b2):
    """The canonical pair (rule, b1, b2) of `strand` of `kind`: MDATK's _pair_fit; for the triplex
    kinds the tables' pose of their canonical triples, or the H-bond fit on the kind's nucleosides."""
    from mdatk import builder as B

    if kind not in EXP_KINDS:
        return B._pair_fit(rule, b1, b2)
    if (rule, b1, b2) in TABLE_PAIRS:
        return table_pose(rule, b1, b2)
    return _cached(("pair", kind, strand, rule, b1, b2),
                   lambda: hbond_pair_fit(rule, b1, b2, kind_core(kind, 0, b1), kind_core(kind, strand, b2)))


def _contact_pairs_of(strands, placed):
    """MDATK's _contact_pairs over levels 0 .. CONTACT_LEVELS (MDATK: 0 .. 2)."""
    import numpy as np
    from mdatk import builder as B

    q = [(0.0, 58.0, 174.0)] * len(strands)
    names, owner, _ = B._level_atoms(placed, strands, q, B._Screw((0, 0, 0), (0, 0, 1), 30.0, 3.0))
    bases = [strands[s].base for s in owner]
    cls = np.array([B._atom_class(b, a) for a, b in zip(names, bases)])
    hbonded = set()
    for s, st in enumerate(strands[1:], start=1):
        for x, y, *_ in B._PAIR_HBONDS.get((st.rule, strands[0].base, st.base), ()):
            hbonded |= {((0, x), (s, y)), ((s, y), (0, x))}
    rows = []
    for i, (a, s) in enumerate(zip(names, owner)):
        for j, (b, t) in enumerate(zip(names, owner)):
            for level in range(CONTACT_LEVELS + 1):
                if level == 0 and j <= i:
                    continue
                if s == t and level == 0:
                    path = B._bond_path(bases[i], a, b)
                elif s == t and level == 1:
                    path = (B._link_path(bases[i], a, bases[j], b) if strands[s].step > 0
                            else B._link_path(bases[j], b, bases[i], a))
                else:
                    path = 99
                if path < 4 or (level == 0 and ((s, a), (t, b)) in hbonded):
                    continue
                limit = B._FIT_LIMITS[cls[i], cls[j]]
                if limit > 0:
                    rows.append((i, j, level, limit, cls[i] < 4 and cls[j] < 4, cls[i] in (4, 5) and cls[j] in (4, 5)))
    i, j, level, limit, heavy, hh = (np.array(c) for c in zip(*rows))
    return i.astype(int), j.astype(int), level.astype(int), limit.astype(float), heavy.astype(bool), hh.astype(bool)


def _distances(pairs, placed, strands, q, screw):
    import numpy as np
    from mdatk import builder as B

    i, j, level, *_rest = pairs
    xyz = B._level_atoms(placed, strands, q, screw)[2]
    levels = np.stack([xyz] + [screw.apply(xyz, k) for k in range(1, CONTACT_LEVELS + 1)])
    return np.linalg.norm(xyz[i] - levels[level, j], axis=-1)


def _perp_basis(a):
    import numpy as np

    u1 = np.cross(a, [1.0, 0.0, 0.0])
    u1 = u1 / np.linalg.norm(u1)
    return u1, np.cross(a, u1)


def exp_screw(kind, p):
    """The screw of parameters p[:6] about the tables' screw of `kind`: axis offset (Å, two components
    in the plane through its origin), tilts (degrees), twist and rise changes (gen_chem_tables
    step_parameters' definitions)."""
    import numpy as np
    from mdatk import builder as B

    sc = triplex()["screws"][kind]
    a0, o0 = np.array(sc["axis"], dtype=np.float64), np.array(sc["origin"], dtype=np.float64)
    u1, u2 = _perp_basis(a0)
    t1, t2 = math.radians(p[2]), math.radians(p[3])
    axis = math.cos(t1) * math.cos(t2) * a0 + math.sin(t2) * u1 - math.sin(t1) * math.cos(t2) * u2
    origin = o0 + p[0] * u1 + p[1] * u2
    return B._Screw(tuple(origin), tuple(axis), float(sc["twist"] + p[4]), float(sc["rise"] + p[5]))


def _exp_level(kind, p):
    """The reference level of an experimental kind at parameters p: strand 1's core as it is, each
    partner's core on its tables pose moved by p's rotation (degrees, about its C1') and shift."""
    import numpy as np
    from mdatk import builder as B

    spec = B._HELICES[kind]
    out = [kind_core(kind, 0, spec[0].base)]
    for s in range(1, len(spec)):
        rot, shift = table_pose(spec[s].rule, spec[0].base, spec[s].base)
        core = B._moved(kind_core(kind, s, spec[s].base), rot, shift)
        v = p[6 + 6 * (s - 1):12 + 6 * (s - 1)]
        r = B._rotvec(np.radians(v[:3]))
        c = core["C1'"]
        out.append({n: (x - c) @ r.T + c + v[3:] for n, x in core.items()})
    return out


def exp_residuals(kind, p, active=None):
    """The restrained reference fit of an experimental kind: MDATK's link residuals and torsion priors
    of every strand, the screw's parameters over the tables' step spreads, the partner poses' moves
    over the tables' pose spreads, the reference pairs' H-bonds about the tables' means, the
    separation of consecutive base planes about the tables' mean, and the contacts (CONTACT_FIT)
    of a level with itself and the CONTACT_LEVELS levels after it (active: the indices of the contact
    pairs scored; None: all)."""
    import numpy as np
    from mdatk import builder as B

    spec = B._HELICES[kind]
    ns = len(spec)
    npose = 6 * (ns - 1)
    tab = triplex()
    sc = tab["screws"][kind]
    screw = exp_screw(kind, p)
    placed = _exp_level(kind, p)
    q = [p[6 + npose + 3 * s:9 + npose + 3 * s] for s in range(ns)]
    out = []
    for core, st, qs in zip(placed, spec, q):
        here = B._nucleotides(*B._core_arrays([core], st.base), *[np.array([v]) for v in qs])
        names = list(here)
        moved = screw.apply(np.concatenate([here[k] for k in names]), st.step)
        there = {k: moved[i:i + 1] for i, k in enumerate(names)}
        out += [B._link_residuals(here, there), B._torsion_priors(*[np.array([v]) for v in qs])]
    out.append(np.asarray(p[:6], dtype=np.float64) / np.array(sc["sd"]))
    spread = {s: (rot, shift) for s, rot, shift in tab["pose_sd"]}
    for s in range(1, ns):
        v = np.asarray(p[6 + 6 * (s - 1):12 + 6 * (s - 1)], dtype=np.float64)
        rot, shift = spread[EXP_KINDS[kind][s]]
        out.append(np.concatenate([v[:3] / rot, v[3:] / shift]))
    hb = {(rule, x, y): (m, sd) for rule, x, y, m, sd in tab["hbonds"]}
    for s in range(1, ns):
        for x, y, *_ in B._PAIR_HBONDS[(spec[s].rule, spec[0].base, spec[s].base)]:
            m, sd = hb[(spec[s].rule, x, y)]
            out.append(np.array([(np.linalg.norm(placed[0][x] - placed[s][y]) - m) / sd]))
    base_atoms = np.array([c[n] for c, st in zip(placed, spec) for n in B.BASE_NAMES[st.base]])
    out.append(np.array([(plane_separation(base_atoms, screw.apply(base_atoms, 1)) - sc["stack"]) / sc["stack_sd"]]))
    pairs = _cached(("exp pairs", kind), lambda: _contact_pairs_of(spec, _exp_level(kind, np.zeros(6 + npose))))
    if active is not None:
        pairs = tuple(c[active] for c in pairs)
    d = _distances(pairs, placed, spec, q, screw)
    out.append(np.clip(pairs[3] - d, 0.0, None) / B._FIT_CONTACT_WEIGHT)
    return np.concatenate(out)


def _exp_near(kind, p):
    """The indices of the contact pairs within REFERENCE_CUTOFF at p."""
    import numpy as np
    from mdatk import builder as B

    spec = B._HELICES[kind]
    ns = len(spec)
    npose = 6 * (ns - 1)
    pairs = _cached(("exp pairs", kind), lambda: _contact_pairs_of(spec, _exp_level(kind, np.zeros(6 + npose))))
    q = [p[6 + npose + 3 * s:9 + npose + 3 * s] for s in range(ns)]
    d = _distances(pairs, _exp_level(kind, p), spec, q, exp_screw(kind, p))
    return np.nonzero(d < REFERENCE_CUTOFF)[0]


def plane_separation(a, b):
    """Separation of the least-squares planes of two point sets along their mean normal."""
    import numpy as np

    ca, cb = a.mean(axis=0), b.mean(axis=0)
    na, nb = np.linalg.svd(a - ca)[2][2], np.linalg.svd(b - cb)[2][2]
    if na @ nb < 0:
        nb = -nb
    n = (na + nb) / np.linalg.norm(na + nb)
    return abs(float((cb - ca) @ n))


def exp_fit(kind):
    """The restrained reference fit (exp_residuals) of an experimental kind from the REFERENCE_STARTS
    (the lowest cost wins): {"x", "cost", "nfev", "screw", "torsions", "poses" (strand -> (rotation,
    shift) of its reference pair), "contact", "contact_h"}."""
    import numpy as np
    from mdatk import builder as B

    def make():
        spec = B._HELICES[kind]
        ns = len(spec)
        npose = 6 * (ns - 1)
        scale = np.array(list(REFERENCE_SCALE_SCREW) + list(REFERENCE_SCALE_POSE) * (ns - 1) + [10.0] * 3 * ns)
        best = None
        for start in REFERENCE_STARTS:
            # the contact pairs within REFERENCE_CUTOFF scored, those within it at a solution added
            # (in the pairs' order) and the solve repeated from there while the set grows
            x = np.concatenate([np.zeros(6 + npose), np.tile(start, ns)])
            active = _exp_near(kind, x)
            nfev = 0
            for _round in range(REFERENCE_ROUNDS):
                sol = SOLVER["least_squares"](lambda p, a=active: exp_residuals(kind, p, a), x, x_scale=scale,
                                              ftol=SMALL_FIT_TOL, xtol=SMALL_FIT_TOL, gtol=SMALL_FIT_TOL,
                                              max_nfev=SMALL_FIT_MAX_NFEV)
                nfev += sol.nfev
                x = sol.x
                grown = np.union1d(active, _exp_near(kind, x))
                if len(grown) == len(active):
                    break
                active = grown
            sol.nfev = nfev
            if best is None or sol.cost < best.cost - 1e-9:
                best = sol
        x = best.x
        screw = exp_screw(kind, x)
        placed = _exp_level(kind, x)
        poses = {}
        for s in range(1, ns):
            rot0, shift0 = table_pose(spec[s].rule, spec[0].base, spec[s].base)
            v = x[6 + 6 * (s - 1):12 + 6 * (s - 1)]
            r = B._rotvec(np.radians(v[:3]))
            c = shift0 + rot0 @ kind_core(kind, s, spec[s].base)["C1'"]
            poses[s] = (r @ rot0, (shift0 - c) @ r.T + c + v[3:])
        torsions = tuple(tuple(float(v) for v in x[6 + npose + 3 * s:9 + npose + 3 * s]) for s in range(ns))
        pairs = _cached(("exp pairs", kind), lambda: _contact_pairs_of(spec, _exp_level(kind, np.zeros(6 + npose))))
        d = _distances(pairs, placed, spec, [np.array(t) for t in torsions], screw)
        return {"x": x, "cost": float(best.cost), "nfev": int(best.nfev), "status": int(best.status),
                "screw": screw, "torsions": torsions, "poses": poses,
                "contact": float(d[pairs[4]].min()), "contact_h": float(d[pairs[5]].min())}
    return _cached(("exp fit", kind), make)


def reference_pose(kind, strand):
    """(rotation, shift) of the reference partner of `strand` (its rule's pair with the strand-1
    reference base): MDATK's _pair_fit, for an experimental kind its fitted pose."""
    from mdatk import builder as B

    spec = B._HELICES[kind]
    if kind in EXP_KINDS:
        return exp_fit(kind)["poses"][strand]
    return B._pair_fit(spec[strand].rule, spec[0].base, spec[strand].base)


def _patch_ch18():
    import numpy as np
    from mdatk import builder as B

    # contacts over the whole helix (CONTACT_LEVELS levels after a level; MDATK: two)
    @functools_cache
    def _contact_pairs(kind):
        return _contact_pairs_of(B._HELICES[kind], B._reference_level(kind)[0])

    def _level_distances(kind, placed, strands, q, screw):
        pairs = B._contact_pairs(kind)
        return _distances(pairs, placed, strands, q, screw), pairs

    B._contact_pairs = _contact_pairs
    B._level_distances = _level_distances

    helix_fit = B._helix_fit

    def _helix_fit(kind):
        if kind in EXP_KINDS:
            f = exp_fit(kind)
            return B._HelixFit(f["screw"], f["torsions"], f["cost"], f["contact"], f["contact_h"])
        return helix_fit(kind)

    B._helix_fit = _helix_fit

    reference_level = B._reference_level

    def _reference_level(kind):
        if kind not in EXP_KINDS:
            return reference_level(kind)
        spec = B._HELICES[kind]
        placed = [kind_core(kind, 0, spec[0].base)] + [B._moved(kind_core(kind, s, spec[s].base),
                                                               *reference_pose(kind, s)) for s in range(1, len(spec))]
        t1 = placed[0]
        n = -B._base_normal(t1, spec[0].base)
        e1 = B._unit(np.cross(n, [1.0, 0.0, 0.0]))
        centre = np.mean([p["C1'"] for p in placed], axis=0)
        return placed, centre, n, e1, np.cross(n, e1)

    B._reference_level = _reference_level

    level_pose = B._level_pose

    def _level_pose(kind, strand, level):
        if kind not in EXP_KINDS:
            return level_pose(kind, strand, level)
        return _cached(("level pose", kind, strand, level), lambda: _exp_level_pose(kind, strand, level))

    def _exp_level_pose(kind, strand, level):  # MDATK's _level_pose on the kind's nucleosides
        from scipy.optimize import least_squares  # noqa: F401  (SOLVER)

        spec = B._HELICES[kind]
        rule, a, b = spec[strand].rule, level[0], level[strand]
        if B._canonical(rule, a, b):
            if (a, b) == (spec[0].base, spec[strand].base):
                return reference_pose(kind, strand)
            return kind_pair_fit(kind, strand, rule, a, b)
        rot0, shift0 = reference_pose(kind, strand)
        others = [kind_core(kind, 0, a)] + [B._moved(kind_core(kind, t, level[t]), *_level_pose(kind, t, level[:t + 1]))
                                            for t in range(1, strand)]
        keep = set(B._RIGID) | {"C1'"}

        def atoms(coords, base):
            names = [x for x in coords if x in keep or (x in B._base_slots(base) and x != "H3+")]
            return np.array([coords[x] for x in names]), np.array([B._atom_class(base, x) for x in names])

        fixed = [atoms(c, base) for c, base in zip(others, level)]
        fixed_xyz = np.vstack([f[0] for f in fixed])
        fixed_cls = np.concatenate([f[1] for f in fixed])
        start = B._moved(kind_core(kind, strand, b), rot0, shift0)
        mobile, mobile_cls = atoms(start, b)
        limit = B._FIT_LIMITS[mobile_cls[:, None], fixed_cls[None, :]]
        planar = np.array([start[x] for x in B.BASE_NAMES[b]])
        ring1, n1 = B._base_plane(others[0], a)
        centre = start["C1'"]

        def move(p, x):
            return (x - centre) @ B._rotvec(np.radians(p[:3])).T + centre + p[3:]

        def residuals(p):
            d = np.linalg.norm(move(p, mobile)[:, None] - fixed_xyz[None], axis=-1)
            return np.concatenate([p[:3] / 10.0, p[3:] / 0.5, ((move(p, planar) - ring1) @ n1) / 0.25,
                                   (np.clip(limit - d, 0.0, None) / 0.05).ravel()])

        best = None
        for turn in (0.0, 25.0, -25.0, 50.0, -50.0):
            sol = SOLVER["least_squares"](residuals, np.concatenate([np.degrees(n1 * math.radians(turn)), np.zeros(3)]),
                                          ftol=SMALL_FIT_TOL, xtol=SMALL_FIT_TOL, gtol=SMALL_FIT_TOL,
                                          max_nfev=SMALL_FIT_MAX_NFEV)
            if best is None or sol.cost < best.cost - 1e-9:
                best = sol
        r = B._rotvec(np.radians(best.x[:3]))
        return r @ rot0, (shift0 - centre) @ r.T + centre + best.x[3:]

    B._level_pose = _level_pose

    helix_strands = B._helix_strands

    def _helix_strands(kind, sequences):
        if kind not in EXP_KINDS:
            return helix_strands(kind, sequences)
        CLOSURE_CONTEXT["kind"] = kind
        fit = B._helix_fit(kind)
        spec = B._HELICES[kind]
        screw = fit.screw
        seq = sequences[0]
        n = len(seq)
        placed = [[{x: screw.apply(v, k) for x, v in kind_core(kind, 0, a).items()} for k, a in enumerate(seq)]]
        level_bases = [seq]
        at_level = [(a,) for a in seq]
        levels, uniform = [], True
        for s, (st, other) in enumerate(zip(spec[1:], sequences[1:]), start=1):
            strand, lvls = [], []
            for k, a in enumerate(seq):
                b = other[k] if st.step > 0 else other[n - 1 - k]
                at_level[k] += (b,)
                canonical = B._canonical(st.rule, a, b)
                rot, shift = B._level_pose(kind, s, at_level[k])
                uniform = uniform and (a, b) == (spec[0].base, st.base)
                strand.append({x: screw.apply(v, k) for x, v in B._moved(kind_core(kind, s, b), rot, shift).items()})
                lvls.append(B._Level(k, a + b, canonical, B._PAIR_HBONDS[(st.rule, a, b)] if canonical else ()))
            placed.append(strand)
            level_bases.append("".join(lv.bases[1] for lv in lvls))
            levels.append(lvls)
        torsions = [np.tile(q, (n, 1)) for q in fit.torsions]
        if not uniform:
            torsions, placed = B._close_levels(placed, level_bases, [st.step for st in spec], torsions)
        strands = []
        for st, nucleosides, bases, q in zip(spec, placed, level_bases, torsions):
            o = slice(None) if st.step > 0 else slice(None, None, -1)
            strands.append(B._strand_heavy(list(nucleosides)[o], bases[o], np.asarray(q)[o]))
        return strands, levels, fit

    B._helix_strands = _helix_strands

    # the refusals, the experimental kinds' REMARK and note, the base-base H-bond note
    check_helix = B.check_helix

    def check_helix_refusing(sequence, partner="auto", pairing="wc", triplex=None, third="auto"):
        try:
            p = B.check_pairing(pairing)
            motif = B.check_triplex(triplex)
            _seq, other, _notes = B.resolve_partner(sequence, partner, p)
            valid = other is not None and not (motif and p != "wc")
        except B.BuildError:
            valid = False
        if valid and refusal(p, motif):
            return [("error", refusal(p, motif))]
        return check_helix(sequence, partner, pairing, triplex, third)

    B.check_helix = check_helix_refusing
    build_helix = B.build_helix

    def build_helix_refusing(sequence, partner="auto", chirality="D", chains=("A", "B"), start=(1, 1), scheme=None,
                             terminals=True, hydrogens=False, phosphate5=False, center=True, notes=None,
                             pairing="wc", triplex=None, third="auto"):
        p = B.check_pairing(pairing)
        motif = B.check_triplex(triplex)
        seq, other, _msgs = B.resolve_partner(sequence, partner, p)
        hand = (chirality or "D").strip().upper()
        if hand in ("D", "L") and other is not None and not (motif and p != "wc") and refusal(p, motif):
            raise B.BuildError(refusal(p, motif))
        mine = []
        text = build_helix(sequence, partner, chirality, chains, start, scheme, terminals, hydrogens, phosphate5,
                           center, mine, pairing, triplex, third)
        kind = (p if p != "wc" else "triplex-" + motif) if (p != "wc" or motif) else None
        if kind in EXP_KINDS:
            text = text.replace("REMARK 951 NUCLEOTIDE OF PDB 4JRD,", "REMARK 951 PDB 3P22 4PLX 6SVS 7JNH,")
            mine = [m.replace("nucleotides (C3'-endo, anti) from PDB 4JRD with backbone torsions fitted to close the "
                              "O3'-P links (reference helix:",
                              "nucleotides (C3'-endo, anti), base pairs and stacking from the RNA triple helices of "
                              f"PDB {EXP_SOURCES} (mean U·A-U triples), the helix, the pairs and the "
                              "backbone torsions fitted to them and to close the O3'-P links (reference helix:")
                    for m in mine]
        if other is not None:
            note = hbond_note(text, seq, other, p, motif, third)
            if note:
                built = next(i for i, m in enumerate(mine) if m.startswith("Built "))
                mine.insert(built, note)
        if notes is not None:
            notes.extend(mine)
        return text

    B.build_helix = build_helix_refusing


def functools_cache(fn):
    import functools

    return functools.lru_cache(maxsize=None)(fn)


CLOSURE_CONTEXT = {}

_DONORS = {"A": ("N6",), "G": ("N1", "N2"), "C": ("N4",), "U": ("N3",)}
_ACCEPTORS = {"A": ("N1", "N3", "N7"), "G": ("O6", "N3", "N7"), "C": ("O2", "N3"), "U": ("O2", "O4")}


def _hbonded(r1, b1, r2, b2):
    import numpy as np

    for donors, acceptors, d_res, a_res in ((_DONORS[b1], _ACCEPTORS[b2], r1, r2),
                                            (_DONORS[b2], _ACCEPTORS[b1], r2, r1)):
        for x in donors:
            for y in acceptors:
                if x in d_res and y in a_res and np.linalg.norm(d_res[x] - a_res[y]) < HBOND_NOTE:
                    return True
    return False


def hbond_note(text, seq, other, pairing, motif, third):
    """The warning of the base pairs (strand-1 base with each partner at its level) whose bases have
    no donor-acceptor pair (N, O) within HBOND_NOTE, or None."""
    import numpy as np
    from mdatk import builder as B
    from mdatk import pdbtools as pt

    s = pt.read_pdb(text)
    segs = s.segments
    n = len(seq)
    rows = []
    partners = [(1, B.PARALLEL[pairing])]
    if motif:
        partners.append((2, motif == "hoogsteen"))
    for k in range(n):
        r1 = segs[0][k]
        b1 = seq[k]
        c1 = {a.name: np.array(a.xyz) for a in r1.atoms}
        for idx, par in partners:
            j = k if par else n - 1 - k
            r2 = segs[idx][j]
            b2 = pt.parse_residue_name(r2.name).base
            c2 = {a.name: np.array(a.xyz) for a in r2.atoms}
            if not _hbonded(c1, b1, c2, b2):
                rows.append(f"{r1.chain}:{r1.resseq} {b1} · {r2.chain}:{r2.resseq} {b2}")
    if not rows:
        return None
    return (f"WARNING: {len(rows)} base pair(s) without a base-base H-bond (no N/O donor-acceptor pair within "
            f"{HBOND_NOTE:.1f} Å): " + ", ".join(rows[:8]) + (" ..." if len(rows) > 8 else ""))


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
