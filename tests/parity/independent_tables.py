"""MDAnalysis with MDBridge's independent data tables patched in.

    import independent_tables
    independent_tables.apply()        # before the first Universe / selection / PDBWriter

    python3 -I -B tests/parity/independent_tables.py   # check the CCD reference names (below)

MDBridge's residue-name, atom-name, element and mass tables are rebuilt from public sources
(src/core/include/mdbridge/core/detail/data_tables_independent.inc, written by
tools/gen_data_tables.py from data/sources/; see data/sources/ATTRIBUTION.md) instead of being copied from MDAnalysis. The generators whose
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
  atomelements (the name-only force-field atom name -> element table, keyed by the upper-case
  name; virtual sites map to "", whose mass is unknown, i.e. 0.0); elements (the guessable one-
  and two-letter elements, upper case as the guesser compares them). The generator's collision
  rule has already removed the entries that would make CCD standard-residue names metals;
* PDB writer (MDAnalysis.coordinates.PDB.PDBWriter): ions (monatomic CCD ion names, after the
  collision rule), special_hg,
  special_cl (CCD components whose HG / CL atom is mercury / chlorine), include_pairs and
  exclude_pairs empty (no independent counterpart).

On top of that, residue_rule_universe() and residue_rule_pdb_text() apply MDBridge's own
residue-scoped ion rule, which MDAnalysis does not have (see the section below).

Not patched (no independent counterpart, not used by MDBridge): vdwradii (bond guessing),
SYBYL2SYMB (MOL2), the raw TABLE_* strings (only read when MDAnalysis is imported).

Run as a script, it checks the result with the patched MDAnalysis: it guesses the element of every
reference name of the collision rule (gen_data_tables.reference_names(): the CCD standard-residue
atom names and their v2.3 spellings, ...) and lists those whose guess is not their CCD element,
and those that the PDB writer would start in the wrong column (column 13 for two-letter elements
and names of four characters, column 14 otherwise). Exit status 1 if any is found.
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
    one, _ = g.guessable_elements()
    two, _ = g.guessable_two_letter()
    names, _ = g.name_elements()
    ions, _ = g.writer_ion_names()
    hg, cl = g.hg_cl_residues()
    elements = g.elements()
    atomelements = {name: element.upper() for name, element, _origins in names}
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
        "ions": tuple(ions),
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


# --- MDBridge's residue-scoped ion rule (not MDAnalysis behaviour) ---------------------------------
#
# MDBridge improves on MDAnalysis' name-only guessing with gen_data_tables.residue_atom_elements()
# (kIndResidueAtomElements; semantics in that generator's docstring and tests/data/README.md).
# MDAnalysis cannot express it, so the generators post-process patched MDAnalysis' results:
#
# * residue_rule_universe(u, path): for every atom whose type MDAnalysis guessed from its name
#   (the topology file gave none, or an empty one) and whose (resname, name) in upper case is in
#   the table, the type becomes the element symbol in upper case; where MDAnalysis took the mass
#   from the type (no elements attribute, no masses in the file), the mass becomes that element's
#   (the patched mass table). Everything else is left as MDAnalysis made it.
# * residue_rule_pdb_text(raw): in a file PDBWriter wrote, an ATOM / HETATM name field written from
#   column 14 (" " + name, name of 2 or 3 characters) whose (resname, name) in upper case is in the
#   table with a two-letter element is rewritten from column 13. The residue name is read from
#   columns 18-21, where PDBWriter writes the first four characters of it.

_residue_elements = None


def residue_elements() -> dict[tuple[str, str], str]:
    global _residue_elements
    if _residue_elements is None:
        _residue_elements = {(res, atom): el for res, atom, el, _ in load_generator().residue_atom_elements()}
    return _residue_elements


def residue_rule_universe(u, path, **kwargs) -> None:
    """Apply the residue rule to the guessed types / masses of u, built from path (kwargs as given
    to the Universe)."""
    import MDAnalysis as mda
    import numpy as np
    from MDAnalysis.exceptions import NoDataError
    from MDAnalysis.guesser import tables as T

    apply()
    table = residue_elements()
    raw = mda.Universe(path, to_guess=(), **kwargs)  # what the file itself gives
    n = len(u.atoms)
    try:
        given = raw._topology.types
        guessed = np.asarray(given.are_values_missing(given.values), dtype=bool)
    except (AttributeError, NoDataError):
        guessed = np.ones(n, dtype=bool)
    mass_from_type = not _has(raw.atoms, "masses") and not _has(u.atoms, "elements")
    try:
        types = np.array(u.atoms.types, dtype=object)
    except NoDataError:
        return
    masses = np.array(u.atoms.masses, dtype=np.float64) if _has(u.atoms, "masses") else None
    changed = False
    for i, atom in enumerate(u.atoms):
        if not guessed[i]:
            continue
        element = table.get((str(atom.resname).upper(), str(atom.name).upper()))
        if element is None:
            continue
        types[i] = element.upper()
        if mass_from_type and masses is not None:
            masses[i] = T.masses[element]
        changed = True
    if changed:
        u.atoms.types = types
        if masses is not None:
            u.atoms.masses = masses


def _has(group, attribute: str) -> bool:
    from MDAnalysis.exceptions import NoDataError
    try:
        getattr(group, attribute)
        return True
    except (AttributeError, NoDataError):
        return False


def residue_rule_pdb_text(raw: bytes) -> bytes:
    """Apply the residue rule to the atom-name columns of a PDB file written by PDBWriter."""
    table = residue_elements()
    text = raw.decode("utf-8")
    lines = text.split("\n")
    for k, line in enumerate(lines):
        if not line.startswith(("ATOM  ", "HETATM")) or len(line) < 21:
            continue
        field, resname = line[12:16], line[17:21].strip()
        name = field.strip()
        if field != " " + name.ljust(3) or not 2 <= len(name) <= 3:
            continue
        element = table.get((resname.upper(), name.upper()))
        if element is not None and len(element) == 2:
            lines[k] = line[:12] + name.ljust(4) + line[16:]
    return "\n".join(lines).encode("utf-8")


def reference_mismatches() -> tuple[list[tuple[str, str, str]], list[tuple[str, str, str]]]:
    """With the patched MDAnalysis: (name, CCD element, guessed type) of the reference names whose
    guess is not their element, and (name, CCD element, written field) of those the PDB writer
    places in the wrong column."""
    import warnings

    from MDAnalysis.coordinates.PDB import PDBWriter
    from MDAnalysis.guesser.default_guesser import DefaultGuesser

    apply()
    guesser = DefaultGuesser(None)
    writer = PDBWriter.__new__(PDBWriter)  # _deduce_PDB_atom_name reads only class data
    guesses, columns = [], []
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for name, els in sorted(load_generator().reference_names().items()):
            for element in sorted(els):
                guess = guesser.guess_atom_element(name)
                if guess.upper() != element:
                    guesses.append((name, element, guess))
                field = writer._deduce_PDB_atom_name(name, "")
                column13 = len(name) >= 4 or len(element) == 2
                if (field[:1] != " ") != column13:
                    columns.append((name, element, field))
    return guesses, columns


def main() -> int:
    guesses, columns = reference_mismatches()
    for name, element, guess in guesses:
        print(f"guess   {name!r}: CCD element {element}, guessed {guess!r}")
    for name, element, field in columns:
        print(f"column  {name!r}: CCD element {element}, written {field!r}")
    print(f"{len(guesses)} guesses and {len(columns)} writer placements differ from the CCD")
    return 1 if guesses or columns else 0


if __name__ == "__main__":
    import sys
    import warnings

    warnings.filterwarnings("ignore")
    sys.exit(main())
