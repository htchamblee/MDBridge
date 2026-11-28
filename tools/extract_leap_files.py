"""Check the AmberTools file and unit names the Amber step types rely on; write data/sources/amber/leap_files.tsv.

    python3 -I tools/extract_leap_files.py AMBERTOOLS PARMED

AMBERTOOLS is an AmberTools tree with dat/leap (an $AMBERHOME, or the conda-forge package
ambertools-26.0 extracted; data/sources/amber/README.md names the build), PARMED a ParmEd source tree
(the PyPI sdist parmed-4.3.1 extracted). For every name of NAMES below the file or unit is looked
up and the row records where it was found and the SHA-256 of that file. Only names (facts) are
written; no file content is copied. The run fails when a name is missing, or when a name recorded
as absent (an empty value) exists after all.

Checks per water model: the leaprc file exists in dat/leap/cmd and names its water unit WAT and HOH; the solvent box unit is
defined in a library the leaprc loads (`loadOff`), its residues are named WAT (the strip masks'
`:WAT`) and it has extra points (atom type EP) exactly when the model is a four-site one; the 12-6-4
ion frcmod is in dat/leap/parm; the ParmEd `add12_6_4 watermodel` name is one ParmEd accepts
(parmed/tools/actions.py `_supported_wms`). The ions of the strip masks are units of
lib/atomic_ions.lib. Standard library only.
"""

from __future__ import annotations

import glob
import hashlib
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "data", "sources", "amber", "leap_files.tsv")

# (table, id, field, value). An empty value records that AmberTools / ParmEd has none.
FORCE_FIELDS = [("ol3cp", "leaprc.RNA.OL3"), ("ol3", "leaprc.RNA.OL3"), ("ljbb", "leaprc.RNA.LJbb"),
                ("roc", "leaprc.RNA.ROC"), ("yil", "leaprc.RNA.YIL"), ("shaw", "leaprc.RNA.Shaw")]
# id -> (leaprc, box unit, ParmEd add12_6_4 water model, 12-6-4 ion frcmod, four-site)
WATERS = {
    "spce": ("leaprc.water.spce", "SPCBOX", "SPCE", "frcmod.ions234lm_1264_spce", False),
    "tip3p": ("leaprc.water.tip3p", "TIP3PBOX", "TIP3P", "frcmod.ions234lm_1264_tip3p", False),
    "opc": ("leaprc.water.opc", "OPCBOX", "OPC", "frcmod.ionslm_1264_opc", True),
    "opc3": ("leaprc.water.opc3", "OPC3BOX", "OPC3", "frcmod.ionslm_1264_opc3", False),
    "tip4pew": ("leaprc.water.tip4pew", "TIP4PEWBOX", "TIP4PEW", "frcmod.ions234lm_1264_tip4pew", True),
    "fb3": ("leaprc.water.fb3", "FB3BOX", "FB3", "frcmod.ionslm_1264_fb3", False),
    "tip4pd": ("leaprc.water.tip4pd", "TIP4PDBOX", "", "", True),
}
# TIP4P/2005: AmberTools has no files of its own; the tip4p2005 helper writes TIP4P2005BOX from the
# TIP4P-Ew box and the build loads the TIP4P-Ew leaprc
TIP4P2005 = ("leaprc.water.tip4pew", "TIP4P2005BOX", "TIP4PEWBOX")
STRIP_ALIASES = ("WAT", "HOH")  # names every water leaprc gives its water unit (`WAT = SPC` ...)
STRIP_IONS = ("Na+", "Cl-", "K+", "MG")  # the ion units of the default strip / solvent masks
# AmberTools files the helpers and the build script read at run time (dat/leap/...)
HELPER_INPUTS = [("cp-lib", "lib/RNA.lib"), ("cp-lib", "parm/parm10.dat"), ("cp-lib", "parm/lj_1264_pol.dat"),
                 ("tip4p2005", "lib/tip4pewbox.off"), ("tip4p2005", "parm/frcmod.tip4pew")]


def sha256(path):
    with open(path, "rb") as fh:
        return hashlib.sha256(fh.read()).hexdigest()


def off_units(path):
    """{unit: {"types": set of atom types, "residues": set of residue names}} of an OFF library."""
    units, unit, block = {}, None, None
    with open(path, encoding="latin-1") as fh:
        for line in fh:
            if line.startswith("!"):
                m = re.match(r"!entry\.(.+?)\.unit\.(atoms|residues) table", line)
                unit, block = (m.group(1), m.group(2)) if m else (None, None)
                if unit:
                    units.setdefault(unit, {"types": set(), "residues": set()})
                continue
            if unit and block:
                fields = line.split('"')
                if len(fields) >= 3:
                    if block == "atoms" and len(fields) >= 5:
                        units[unit]["types"].add(fields[3])
                    elif block == "residues":
                        units[unit]["residues"].add(fields[1])
    return units


def main(argv=None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if len(argv) != 2:
        raise SystemExit(__doc__.split("\n\n")[1])
    leap = os.path.join(argv[0], "dat", "leap")
    actions = os.path.join(argv[1], "parmed", "tools", "actions.py")
    with open(actions, encoding="utf-8") as fh:
        m = re.search(r"_supported_wms = \(([^)]*)\)", fh.read())
    if not m:
        raise SystemExit(f"{actions}: no _supported_wms")
    parmed_models = re.findall(r"'([^']+)'", m.group(1))
    parmed_rel = "parmed/tools/actions.py"
    parmed_sha = sha256(actions)
    rows = []

    def found(rel):
        path = os.path.join(leap, rel)
        if not os.path.isfile(path):
            raise SystemExit(f"dat/leap/{rel}: not found")
        return path

    def unit_in(lib_names, unit):
        for lib in lib_names:
            units = off_units(found(f"lib/{lib}"))
            if unit in units:
                return lib, units[unit]
        raise SystemExit(f"unit {unit} not in {', '.join(lib_names)}")

    def leaprc_libs(rc):
        with open(found(f"cmd/{rc}"), encoding="latin-1") as fh:
            text = fh.read()
        for alias in STRIP_ALIASES:
            if not re.search(rf"^{alias} = ", text, re.M):
                raise SystemExit(f"{rc}: does not set {alias}")
        return re.findall(r"^loadOff (\S+)", text, re.M)

    for fid, rc in FORCE_FIELDS:
        rows.append(("forcefield", fid, "leaprc", rc, f"dat/leap/cmd/{rc}", sha256(found(f"cmd/{rc}"))))
    for wid, (rc, box, model, frcmod, four) in WATERS.items():
        rows.append(("water", wid, "leaprc", rc, f"dat/leap/cmd/{rc}", sha256(found(f"cmd/{rc}"))))
        libs = leaprc_libs(rc)
        # the library the leaprc loads last that defines the box (later loadOff calls replace units)
        lib, unit = unit_in(list(reversed(libs)), box)
        if unit["residues"] != {"WAT"}:
            raise SystemExit(f"{box}: residues {sorted(unit['residues'])}, expected WAT")
        if ("EP" in unit["types"]) != four:
            raise SystemExit(f"{box}: extra points {'EP' in unit['types']}, expected {four}")
        rows.append(("water", wid, "box", box, f"dat/leap/lib/{lib}", sha256(found(f"lib/{lib}"))))
        rows.append(("water", wid, "four_site", "1" if four else "0", f"dat/leap/lib/{lib}",
                     sha256(found(f"lib/{lib}"))))
        if model:
            if model not in parmed_models:
                raise SystemExit(f"ParmEd has no add12_6_4 water model {model}")
            rows.append(("water", wid, "parmed_1264", model, parmed_rel, parmed_sha))
        else:
            stem = rc.split(".")[-1]
            if stem.upper() in parmed_models:
                raise SystemExit(f"ParmEd has a 12-6-4 water model {stem.upper()} after all")
            rows.append(("water", wid, "parmed_1264", "", parmed_rel, parmed_sha))
        if frcmod:
            rows.append(("water", wid, "frcmod_1264", frcmod, f"dat/leap/parm/{frcmod}",
                         sha256(found(f"parm/{frcmod}"))))
        else:
            stem = rc.split(".")[-1]
            hits = glob.glob(os.path.join(leap, "parm", f"frcmod.*1264*{stem}"))
            if hits:
                raise SystemExit(f"a 12-6-4 frcmod for {stem} exists: {hits}")
            rows.append(("water", wid, "frcmod_1264", "", f"dat/leap/parm/frcmod.*1264*{stem}", ""))
    rc, box, base = TIP4P2005
    rows.append(("water", "tip4p2005", "leaprc", rc, f"dat/leap/cmd/{rc}", sha256(found(f"cmd/{rc}"))))
    lib, unit = unit_in(list(reversed(leaprc_libs(rc))), base)
    for name in glob.glob(os.path.join(leap, "lib", "*")):
        if name.endswith((".lib", ".off")) and box in off_units(name):
            raise SystemExit(f"{box} exists in AmberTools ({name})")
    rows.append(("water", "tip4p2005", "box", box, "written by the tip4p2005 helper from "
                 f"dat/leap/lib/tip4pewbox.off ({base})", ""))
    rows.append(("water", "tip4p2005", "four_site", "1", f"dat/leap/lib/{lib}", sha256(found(f"lib/{lib}"))))
    rows.append(("water", "tip4p2005", "parmed_1264", "", parmed_rel, parmed_sha))
    rows.append(("water", "tip4p2005", "frcmod_1264", "", "dat/leap/parm (none for TIP4P/2005)", ""))
    for alias in STRIP_ALIASES:
        rows.append(("strip", "water", "alias", alias, "dat/leap/cmd/leaprc.water.* (every water model above)", ""))
    for ion in STRIP_IONS:
        lib, _unit = unit_in(["atomic_ions.lib"], ion)
        rows.append(("strip", "ion", "unit", ion, f"dat/leap/lib/{lib}", sha256(found(f"lib/{lib}"))))
    for helper, rel in HELPER_INPUTS:
        rows.append(("helper_input", helper, "file", rel, f"dat/leap/{rel}", sha256(found(rel))))

    with open(OUT, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("table\tid\tfield\tvalue\tfound_in\tsha256\n")
        for row in rows:
            fh.write("\t".join(row) + "\n")
    print(f"{len(rows)} rows")
    return 0


if __name__ == "__main__":
    sys.exit(main())
