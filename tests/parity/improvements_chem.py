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
         MDATK's wrap
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
  D7     dict pubchem asks MDBRIDGE_PUBCHEM_URL (the PUG REST root, default MDATK's PUBCHEM_REST) when
         that variable is set: the local stand-in server of cli_parity.py
  D9(M4) `render` / `preview` fit a shown map without a placement again when auto_fit is on (MDATK's
         ensure_density_fit; in the command line only: the milestone-4 render generators keep D9)
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
    _patch_pdbtools()
    _patch_annotate()
    _patch_builder()
    _patch_density()


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

    # the density fit before rendering (auto_fit) is back: MDATK 0.0.4.1's ensure_density_fit in place
    # of milestone 4's D9 stand-in (improvements_render._d9_no_density_fit, which the milestone-4
    # generators still apply; W2-C retires it there)
    def ensure_density_fit(self, reporter=None):
        from mdatk.progress import Cancelled

        p = self.project
        rep = reporter or self.reporter
        if p.density_active() and len(p.density.transform) != 12:
            if not p.density.auto_fit:
                rep.log("WARNING: the density map has no placement and automatic fitting is off; it is not drawn")
                return
            try:
                self.fit_density(reporter=reporter)
            except Cancelled:
                raise
            except (pipeline.PipelineError, OSError, ValueError) as exc:
                rep.log(f"WARNING: density map not drawn: {exc}")

    pipeline.Pipeline.ensure_density_fit = ensure_density_fit

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
