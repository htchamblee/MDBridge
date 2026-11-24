"""MDAnalysis with MDBridge's independent data tables patched in.

    import independent_tables
    independent_tables.apply()        # before the first Universe / selection / PDBWriter

MDBridge's residue-name, atom-name, element and mass tables are rebuilt from public sources
(src/core/data_tables_independent.inc, written by tools/gen_data_tables.py from data/sources/;
see data/sources/ATTRIBUTION.md) instead of being copied from MDAnalysis. The generators whose
output depends on those tables call apply(), so the committed expected data describes MDAnalysis
2.11.0's behaviour (its algorithms, unchanged) with MDBridge's tables in place of its own.

The tables are loaded with tools/gen_data_tables.py's own functions (from the data/sources/ next
to it), never by parsing the .inc. Only data is replaced, no function. Mutable tables are changed
in place, because modules hold references taken at import time (the parsers' `from
..guesser.tables import SYMB2Z`, LAMMPSParser's `masses as mass_table`, MDAnalysis.topology.tables,
the selection subclasses sharing their parent's residue set); PDBWriter's tuples are replaced on
the class (MultiPDBWriter inherits them). What is patched:

* selection macros (MDAnalysis.core.selection): ProteinSelection.prot_res (protein and backbone),
  NucleicSelection.nucl_res (nucleic, nucleicbackbone, nucleicbase, nucleicsugar),
  WaterSelection.water_res, BackboneSelection.bb_atoms, NucleicBackboneSelection.bb_atoms,
  BaseSelection.base_atoms, NucleicSugarSelection.sug_atoms;
* guessers (MDAnalysis.guesser.tables): Z2SYMB / SYMB2Z (CIAAW symbols, Z 1-118); masses (the
  CIAAW 2024 abridged standard atomic weight, or the mass number of the most stable isotope, of
  every element, under its symbol and its upper-case spelling, as MDAnalysis looks masses up);
  atomelements (the force-field atom name -> element table, keyed by the upper-case name as
  MDAnalysis' guesser looks names up; MDAnalysis has no residue scope, so the residue column of
  the table is dropped; virtual sites map to "", whose mass is unknown, i.e. 0.0); elements (the
  guessable one- and two-letter elements, upper case as the guesser compares them);
* PDB writer (MDAnalysis.coordinates.PDB.PDBWriter): ions (monatomic CCD ion names), special_hg,
  special_cl (CCD components whose HG / CL atom is mercury / chlorine), include_pairs and
  exclude_pairs empty (no independent counterpart).

Not patched (no independent counterpart, not used by MDBridge): vdwradii (bond guessing),
SYBYL2SYMB (MOL2), the raw TABLE_* strings (only read when MDAnalysis is imported).
"""

from __future__ import annotations

import importlib.util
import os

HERE = os.path.dirname(os.path.abspath(__file__))
GENERATOR = os.path.join(os.path.dirname(os.path.dirname(HERE)), "tools", "gen_data_tables.py")

_applied = False


def load_generator():
    """tools/gen_data_tables.py as a module (reads data/sources/ when imported; writes nothing)."""
    spec = importlib.util.spec_from_file_location("gen_data_tables", GENERATOR)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def tables() -> dict:
    """The independent tables, in the shapes MDAnalysis keeps its own."""
    g = load_generator()
    backbone, sugar, base = g.nucleic_sets()
    one, two = g.guessable_elements()
    hg, cl = g.hg_cl_residues()
    elements = g.elements()
    atomelements = {}
    for _residue, name, element, _origin in g.forcefield_name_elements():
        key, value = name.upper(), element.upper()
        if atomelements.setdefault(key, value) != value:
            raise ValueError(f"atom name {name!r} maps to {atomelements[key]!r} and {value!r}")
    masses = {}
    for _z, symbol, _name, mass, _std in elements:
        masses[symbol] = float(mass)
        masses[symbol.upper()] = float(mass)
    return {
        "prot_res": g.protein_residues(),
        "nucl_res": g.nucleic_residues(),
        "water_res": g.water_residues(),
        "backbone": g.protein_backbone(),
        "nucleicbackbone": backbone,
        "nucleicbase": base,
        "nucleicsugar": sugar,
        "z2symb": {z: symbol for z, symbol, _name, _mass, _std in elements},
        "masses": masses,
        "atomelements": atomelements,
        "elements": [e.upper() for e in one + two],
        "ions": tuple(g.monatomic_ion_names()),
        "special_hg": tuple(hg),
        "special_cl": tuple(cl),
    }


def _replace_set(owner, attribute: str, values) -> None:
    current = owner.__dict__[attribute]  # the class's own set, not an inherited one
    if not isinstance(current, set):
        raise TypeError(f"{owner.__name__}.{attribute} is a {type(current).__name__}, not a set")
    current.clear()
    current.update(values)


def _replace_dict(current: dict, values: dict) -> None:
    current.clear()
    current.update(values)


def apply() -> None:
    """Patch the imported MDAnalysis in place (idempotent)."""
    global _applied
    if _applied:
        return
    from MDAnalysis.coordinates.PDB import PDBWriter
    from MDAnalysis.core import selection as S
    from MDAnalysis.guesser import tables as T

    t = tables()
    _replace_set(S.ProteinSelection, "prot_res", t["prot_res"])
    _replace_set(S.NucleicSelection, "nucl_res", t["nucl_res"])
    _replace_set(S.WaterSelection, "water_res", t["water_res"])
    _replace_set(S.BackboneSelection, "bb_atoms", t["backbone"])
    _replace_set(S.NucleicBackboneSelection, "bb_atoms", t["nucleicbackbone"])
    _replace_set(S.BaseSelection, "base_atoms", t["nucleicbase"])
    _replace_set(S.NucleicSugarSelection, "sug_atoms", t["nucleicsugar"])

    _replace_dict(T.Z2SYMB, t["z2symb"])
    _replace_dict(T.SYMB2Z, {v: k for k, v in t["z2symb"].items()})
    _replace_dict(T.masses, t["masses"])
    _replace_dict(T.atomelements, t["atomelements"])
    T.elements[:] = t["elements"]

    for name in ("ions", "special_hg", "special_cl"):
        if not isinstance(PDBWriter.__dict__[name], tuple):
            raise TypeError(f"PDBWriter.{name} is not a tuple")
        setattr(PDBWriter, name, t[name])
    PDBWriter.include_pairs = ()
    PDBWriter.exclude_pairs = ()
    _applied = True
