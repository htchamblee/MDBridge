"""DNA / L-DNA residue names, the CCD's own atom names of 0DA / 0DT and the backbone guard (docs/DIFFERENCES.md
NA-01), as test-side patches of MDATK v0.0.4.1, so that the generators record MDBridge's behaviour.

MDBridge follows MDATK 0.1.x here (mdatk/pdbtools.py parse_residue_name, has_backbone_atoms,
standard_atom_name; secstruct.PairDetector, analysis/stacking.py ring_atoms, analysis/structure.py's
trace flag `backbone`, render/molecule/topology.py, render/residues.py one_letter and the sequence
panel), with the residue and atom-name tables taken from the wwPDB Chemical Component Dictionary
extract data/sources/ccd/l_deoxyribonucleotide_atoms.tsv (as tools/gen_data_tables.py reads it):

* parse_residue_name also reads L-DNA: the CCD's L-deoxyribonucleotides 0DA 0DC 0DG 0DT and MDATK's
  LDA LDC LDG LDT (LDA5 / LDA3 / LDAN);
* a residue is a nucleotide by its name only with a backbone atom (a primed name, '*' or P): the
  base-pair detector and the sequence panel (through the trace's per-row `backbone` flags) leave out
  a lanthanum ion LA, radium RA or a detergent LDA;
* the atoms that the CCD names its own way in 0DA / 0DT (N9A ..., C5M) are read by their standard
  names (N9 ..., C7) by the base-pair detector, the stacking rings and the molecule topology.

apply_mdatk() patches the imported MDATK modules (idempotent); improvements_analyses.apply_mdatk()
calls it. parse_residue_name(name) is the improved function for generators that record it alone.
"""

import csv
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
CCD = os.path.join(ROOT, "data", "sources", "ccd", "l_deoxyribonucleotide_atoms.tsv")

_LDNA = re.compile(r"^LD([ACGT])([53N]?)$")


def _ccd_tables():
    with open(CCD, encoding="utf-8", newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    bases = {r["comp_id"]: r["base"] for r in rows}
    own = {}
    for r in rows:
        if r["atom_id"] != r["standard_atom_id"]:
            own.setdefault(r["comp_id"], {})[r["atom_id"]] = r["standard_atom_id"]
    return bases, own


L_DEOXY, STANDARD_NAMES = _ccd_tables()
_original_parse = None
_applied = False


def has_backbone_atoms(atom_names) -> bool:
    return any("'" in n or "*" in n or n.strip() == "P" for n in map(str, atom_names))


def standard_atom_name(resname: str, name: str) -> str:
    return STANDARD_NAMES.get(str(resname).strip().upper(), {}).get(name, name)


def parse_residue_name(name):
    """MDATK 0.0.4.1's parse_residue_name, then the L-DNA names."""
    global _original_parse
    if _original_parse is None:
        from mdatk import pdbtools
        _original_parse = getattr(pdbtools, "_mdbridge_original_parse", pdbtools.parse_residue_name)
    info = _original_parse(name)
    if info is not None:
        return info
    n = name.strip().upper()
    from mdatk.pdbtools import ResidueName
    if n in L_DEOXY:
        return ResidueName(L_DEOXY[n], "dna", "L", "")
    m = _LDNA.match(n)
    return ResidueName(m.group(1), "dna", "L", m.group(2)) if m else None


def standard_index(resname, index):
    """An atom index by standard names; an atom of the standard name wins over a renamed one."""
    out = {name: row for name, row in index.items() if standard_atom_name(resname, name) == name}
    for name, row in index.items():
        out.setdefault(standard_atom_name(resname, name), row)
    return out


def apply_mdatk():
    global _applied
    if _applied:
        return
    _applied = True
    import numpy as np
    from mdatk import pdbtools, secstruct
    from mdatk.analysis import stacking, structure
    from mdatk.render import residues, seqpanels
    from mdatk.render.molecule import topology

    pdbtools._mdbridge_original_parse = pdbtools.parse_residue_name
    secstruct.parse_residue_name = parse_residue_name  # base_letter, canonical_name

    # base pairs: standard atom names, and the backbone guard
    init = secstruct.PairDetector.__init__

    def detector_init(self, resnames, atom_index, cutoff=3.5, wobble=True):
        index = [standard_index(r, idx) for r, idx in zip(resnames, atom_index)]
        # a residue without a backbone atom gets no base (base_letter("") is None)
        names = [r if has_backbone_atoms(idx) else "" for r, idx in zip(resnames, index)]
        init(self, names, index, cutoff=cutoff, wobble=wobble)

    secstruct.PairDetector.__init__ = detector_init

    # stacking rings by standard atom names
    def ring_atoms(residue):
        names = {standard_atom_name(str(residue.resname), str(a.name)): a.ix for a in residue.atoms}
        for ring in (stacking.PURINE, stacking.PYRIMIDINE):
            if all(n in names for n in ring):
                return [names[n] for n in ring]
        ring = stacking.AROMATIC_RESIDUES.get(residue.resname.upper())
        if ring and all(n in names for n in ring):
            return [names[n] for n in ring]
        return None

    stacking.ring_atoms = ring_atoms

    # the trace's backbone flags
    finalize = structure.TraceAnalysis.finalize

    def trace_finalize(self):
        res = finalize(self)
        res.arrays["backbone"] = np.array([has_backbone_atoms(str(n) for n in r.atoms.names)
                                           for r in self.residues], dtype=bool)
        return res

    structure.TraceAnalysis.finalize = trace_finalize

    # the molecule topology reads standard atom names
    topology_init = topology.Topology.__init__

    def topology_init_standard(self, names, elements, resnames, *args, **kwargs):
        names = [standard_atom_name(str(rn), str(n)) for n, rn in zip(names, resnames)]
        topology_init(self, names, elements, resnames, *args, **kwargs)

    topology.Topology.__init__ = topology_init_standard

    # one_letter with the backbone guard; the sequence panel passes the trace's flags
    amino = residues._AMINO

    def one_letter(resname, backbone=None):
        info = parse_residue_name(resname) if backbone is not False else None
        if info is not None:
            return info.base
        return amino.get(resname.strip().upper(), "X")

    residues.one_letter = one_letter
    seqpanels.one_letter = one_letter
    original_residues = seqpanels._ResiduePanel._residues

    def _residues(self):
        original_residues(self)
        flags = self.tr.arrays.get("backbone")
        if flags is not None:
            flags = np.asarray(flags, dtype=bool)
            names = seqpanels.resnames(self.tr)
            self.letters = [one_letter(names[r], bool(flags[r])) for r in self.rows]

    seqpanels._ResiduePanel._residues = _residues
