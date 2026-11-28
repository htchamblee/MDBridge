"""MDBridge's deliberate differences in MDATK v0.0.4.1's workflow layer (docs/specs/m7-plan.md
section 5, WG-01 ... WG-21, and the owner decisions D3, D6, D9, D11, D13, D14 of section 7), as
test-side patches of MDATK, so that make_workflow_cases.py, make_workflow_goldens.py and
cli_parity.py record MDBridge's behaviour. Where a later MDATK fixed the item, its design is followed
and named (MDATK's own code, ported); every other byte stays as MDATK v0.0.4.1 writes it.

apply_mdatk() patches the workflow library (idempotent; it applies improvements_chem.apply_mdatk()
and improvements_citations.apply_mdatk() first: the generators read input files with the chem
layer and cite with MDBridge's table); apply_cli() also the `workflow` commands (the child of
cli_parity.py calls both). mask() is the masking of generated texts (D11). Each patch with its id:

  D11    the program in generated files: "written by MDBridge <VERSION>", "MDBridge only writes this
         script", "# --- settings written by MDBridge ---" (every text MDATK builds from PROJECT_ABBREV
         and __version__ in mdatk/workflows, and embed's settings line); the version is written as
         <VERSION> and every date (YYYY-MM-DD) as <DATE> (mask(); MDATK 0.1.2 tests/golden.py
         normalises the same way, VERSION_RE / DATE_RE). Literal "MDATK" texts (step labels such as
         "MDATK analysis", "MDATK table", messages) stay as MDATK has them
  D13    (WG-11) a float setting that is NaN or infinite is refused by Option.coerce: "buffer: expected
         a finite number, got nan" (inf, -inf), so normalise() falls back to the default with that
         warning and the command line reports it; MDATK 0.1.2 has not changed it
  WG-12  an integer setting given as inf / 1e400: "cannot convert float infinity to integer" (Python's
         OverflowError text, raised as the ValueError normalise and the command line report; MDBridge's
         core::Option::coerce already words it so, without the key, as for NaN)
  WG-01  (MDATK 0.1.0) every hand-picked input path in a shell line of a run / run-separate script
         goes through shell_quote (unchanged text for ordinary paths): words of command lines outside
         comments and here-documents; a DISANG file (connected or chosen) whose path has white space is
         an error of the MD step (the mdin namelist cannot hold it): "Restraints (DISANG): pmemd reads
         DISANG= up to the first space, and this path has one: <path>; give the file a path without
         spaces" (keys stages.<k>.distres, input:restraints)
  WG-02  the master script does not call a step that reads (through connections, directly or not)
         from a step with errors: "#   (not run: reads from <label>, which has errors)" (several:
         "<a>, <b>, which have errors") instead of its bash / run-separately lines; `workflow run`
         (apply_cli) reports each selected step with errors ("error: <label> has errors and was not
         run: <errors>"), does not start the steps that read from it ("<label>: not started (it reads
         from a step with errors)"), runs the others and exits 1 (MDATK 0.1.0 cli._run_plan_scripts)
  WG-03  (MDATK 0.1.2 Plan._check_shared_outputs) an output file that another step writes too is an
         error of the later step: "output 'best': <path> is also written by <label>; give one of them
         another file name" (pairs with the same step number and tag keep MDATK's own error)
  WG-04  (MDATK 0.1.0) a stage in segments with a temperature ramp ramps in its first segment only:
         segments 2-N read <stage>-cont.in, the stage without the ramp (temp0, no TEMP0 schedule)
  WG-05  (MDATK 0.1.0) a stage named like another stage's segment run (<stage>-001 ...) or its
         continuation input (<stage>-cont) is an error naming both: "stages eq (segment 1) and eq-001
         would write the same files (<project>_<NN>_eq-001_…): rename one of them"
  WG-06  (MDATK 0.1.0 count_solute_residues) the solute residue count behind the default restraint
         mask :1-N counts residues as tleap reads them: a new residue wherever the name, chain or number
         changes or a TER intervenes; water names of the PDB preparation (prep_pdb.WATERS) are solvent
  WG-07  (MDATK 0.1.0) FARFAR2 and SimRNA run scripts end with "ERROR: every rna_denovo / SimRNA run
         failed (see above)" and status 1 when every run failed, and warn "WARNING: N of M ... runs
         failed (see above): their models / trajectories are missing" otherwise (a failed run returns 1
         and is counted in a temporary folder)
  WG-08  (MDATK 0.1.0) SimRNA clustering: `REF=$(ls ... | head -n 1 || true)`: a SIGPIPE of ls under
         pipefail no longer ends the script without a message
  WG-09  (MDATK 0.1.0) each script removes the outputs of an earlier run first: the SimRNA clusters and
         models ("rm -f <all>_thrs*_clust* <project>_<NN>_model-[0-9]*_<stamp>.pdb"), the FARFAR2
         models and their list, the RNAMake design PDBs, the VMD movie frames
  WG-10  (MDATK 0.1.0) a secondary-structure file with the dot-bracket line only and no sequence typed:
         "the secondary-structure file has no sequence line: type the sequence" (keys sequence,
         input:<port>) in every step that reads one (amber-restraints, simrna, phenix-rsr)
  WG-19  (MDATK 0.1.0) the FARFAR2 models script annotates every extracted model with DSSR: the loop
         globs <project>_<NN>_model-[0-9]*_<stamp>.pdb (0.0.4.1's model-[0-9][0-9] left out models 100 on)
  WG-20  (MDATK 0.1.0) a SimRNA step without clustering writes its output 'trajectory' itself: its run
         script ends with "# all replica trajectories of all runs", cat runs/r*/*.trafl > <all>.trafl and
         the check that it was written (0.0.4.1 promised the file, and only the clustering script wrote it)
  WG-21  (MDATK 0.1.0 gmx_posre.is_nucleic) the GROMACS position-restraint helper takes the L-DNA residue
         names LDA LDC LDG LDT (and their terminal / N variants) for nucleic acid: an L-DNA strand got no
         restraints; every name recognised before is recognised as before (the helper's only change, D4)
  WG-18  (D6) the structure builder's script starts MDBridge, "$MDBRIDGE" build ..., unless the step's
         Python setting names another interpreter (then MDATK's "<python> -m mdatk" stays); the plan
         summary records the program under each step that runs it (see D3)
  D3     the amber-build L-RNA naming and the restraint dictionaries' renaming run MDBridge instead of
         MDATK's embedded copies of pdbtools.py / dictionaries.py (no such helper file is written):
         "$MDBRIDGE" pdb clean ... / "$MDBRIDGE" dict finish ... with MDATK's arguments, their checks
         "ERROR: mdbridge pdb clean did not write ..." / "mdbridge dict finish"; every script that runs
         MDBridge sets, after its header (and environment lines):
             # MDBridge: the mdbridge program on PATH, unless MDBRIDGE names another
             MDBRIDGE="${MDBRIDGE:-mdbridge}"
         (D6: an AppImage writes its own path instead of mdbridge; not in these cases), and the plan
         summary adds, after the step's "Type:" line (which gets Markdown's two trailing spaces):
             MDBridge: `mdbridge` from PATH (MDBRIDGE in the scripts overrides it)
  D9     (MDATK 0.0.5 / 0.0.5.1 / 0.0.11, physics fixes) amber-build:
         * "Separate HO and O2 Lennard-Jones types first" (split_types) is on by default, and with it
           off a solvated NBfix19 build warns (0.0.5.1's warning instead of 0.0.4.1's info);
         * L-RNA atom types (lrna_types) take TYPE or TYPE:REFERENCE; with microMg / nanoMg every
           L-RNA type gets the Mg2+ pair term of its D-RNA twin (LD and LE as OS, LH as OH by default),
           written after the CP lines ("# non-standard atom types (e.g. an L-RNA library) with the
           terms of their chemical twin", "# LD as OS", changeLJPair ...); a type without a reference is
           an error, the mapping an info; the L-RNA fixes no longer depend on MDATK naming the L
           residues (an input PDB may carry L names already); 0.0.5's messages, except that the
           warning about atom types of extra libraries / parameter files (0.0.5 words it with its
           topology check, which D9 leaves for later) reads: "the <set> pair terms cover the 27
           parmBSC0χOL3 RNA atom types and the L-RNA atom types; other atom types of extra residue
           libraries or parameter files keep the plain combination rule with Mg²⁺", and the CP note
           is 0.0.5's info;
         * the microMg / nanoMg frcmod defines the ion's terms for both type names, Mg2+ and MG;
         * a two-pass build leaves Mg2+ lines (:MG) of the extra ParmEd commands out of the dry pass;
         * the lanthanide trap: with the input PDB present and L residues in it, the names MDATK
           writes for them must be units of the own libraries (Extra residue libraries, loadOff lines
           before loadpdb, units assigned with NAME = ...) or mapped onto one with addPdbResMap;
           otherwise an error with the map line to add (MDATK 0.0.11 _check_lrna_names, its texts with
           "MDATK" as MDATK words them; the warning for a structure of an earlier step is not ported)
  D14    (WG-13) `workflow run` first generates each selected step's files and compares them with
         the files on disk (dates and the version masked); a step whose files differ or are missing
         is reported, "error: step 03: its files on disk were written with other settings; write the
         plan first" (missing: "... its files are not written yet; write the plan first"), and nothing
         runs (exit 1), unless --force (a new option of `workflow run`) runs them anyway

Not patched here: WG-14 (Windows file locks), WG-16 (ending programs a script left running) and
D12 (Windows does not run scripts) are behaviour of the writer and the runner, tested in C++; WG-15
(New stamp) is new_stamp_after() below, used by make_workflow_cases.py; WG-17 (Windows paths) needs a
Windows host to show.
"""

from __future__ import annotations

import datetime as _real_dt
import os
import re
import sys
import types

NAME = "MDBridge"
VERSION = "<VERSION>"
DATE = "<DATE>"
VERSION_RE = re.compile(r"\b(MDBridge|MDATK|mdatk|mdbridge) \d+\.\d+\.\d+(?:\.\d+)?(?:[-+][0-9A-Za-z.]+)?")
DATE_RE = re.compile(r"\b20\d\d-\d\d-\d\d\b")
MDBRIDGE_LINES = ["# MDBridge: the mdbridge program on PATH, unless MDBRIDGE names another",
                  'MDBRIDGE="${MDBRIDGE:-mdbridge}"', ""]
MDBRIDGE_SUMMARY = "MDBridge: `mdbridge` from PATH (MDBRIDGE in the scripts overrides it)"
MDBRIDGE_WORD = '"$MDBRIDGE"'

_applied = False
_applied_cli = False


def mask(text: str) -> str:
    """D11: the version after the program's name and dates, as tests/golden.py of MDATK 0.1.2 masks them."""
    text = VERSION_RE.sub(lambda m: f"{m.group(1)} {VERSION}", text)
    return DATE_RE.sub(DATE, text)


def _sub(text: str, old: str, new: str, count: int = 1) -> str:
    """text with `old` replaced (it must occur at least `count` times: MDATK's text is pinned)."""
    if text.count(old) < count:
        raise AssertionError(f"MDATK's text changed: {old!r} not found")
    return text.replace(old, new)


def apply_mdatk():
    global _applied
    if _applied:
        return
    _applied = True
    here = os.path.dirname(os.path.abspath(__file__))
    if here not in sys.path:
        sys.path.insert(0, here)
    import improvements_chem
    import improvements_citations

    improvements_chem.apply_mdatk()
    from mdatk import citations

    if citations.pipeline_resources is not improvements_citations.pipeline_resources:  # (not applied yet)
        improvements_citations.apply_mdatk()
    _patch_names()
    _patch_coerce()
    _patch_secondary_structure()
    _patch_amber_build()
    _patch_amber_md()
    _patch_modelling()
    _patch_builder()
    _patch_dicts()
    _patch_plan()


# ---------------------------------------------------------------------------------------------
# D11

class _Today:
    def isoformat(self):
        return DATE


class _Date:
    @staticmethod
    def today():
        return _Today()


_DT = types.SimpleNamespace(date=_Date, datetime=_real_dt.datetime, timedelta=_real_dt.timedelta)


def _patch_names():
    from mdatk.workflows import (amber, base, dictstep, farfar2, gromacs, plan, refine, runner,  # noqa: F401
                                 scripts, vmd)

    for module in (amber, base, dictstep, farfar2, gromacs, plan, refine, runner, vmd):
        if hasattr(module, "PROJECT_ABBREV"):
            module.PROJECT_ABBREV = NAME
        if hasattr(module, "__version__"):
            module.__version__ = VERSION
    base._dt = _DT
    plan._dt = _DT
    embed = scripts.embed

    def embed_named(module, constants, call, doc, imports=()):
        return _sub(embed(module, constants, call, doc, imports), "# --- settings written by MDATK ---",
                    f"# --- settings written by {NAME} ---")

    scripts.embed = embed_named
    for module in (amber, gromacs):
        if hasattr(module, "embed"):
            module.embed = embed_named


# ---------------------------------------------------------------------------------------------
# D13 / WG-11, WG-12

def _patch_coerce():
    import math

    from mdatk.workflows import schema

    coerce = schema.Option.coerce

    def coerce_finite(self, value):
        if self.kind in ("int", "float") and value is not None and not isinstance(value, bool) \
                and not (isinstance(value, str) and not value.strip()):
            try:
                number = float(value)
            except (TypeError, ValueError, OverflowError):
                number = None
            if number is not None and not math.isfinite(number):
                if self.kind == "float":  # D13 / WG-11
                    raise ValueError(f"{self.key}: expected a finite number, got {number!r}")
                if math.isinf(number):  # WG-12 (MDATK: an OverflowError nobody catches)
                    raise ValueError("cannot convert float infinity to integer")
        return coerce(self, value)

    schema.Option.coerce = coerce_finite


# ---------------------------------------------------------------------------------------------
# WG-10

def _patch_secondary_structure():
    from mdatk.workflows import amber
    from mdatk.workflows.schema import Issue

    original = amber.secondary_structure

    def secondary_structure(s, ctx, port="secstruct"):
        model, source, issues = original(s, ctx, port)
        value = ctx.value(port)
        if model is not None and not model.chains and value is not None and value.paths \
                and value.data.get("rna_model") is None:
            issues.append(Issue("error", "the secondary-structure file has no sequence line: type the sequence",
                                ("sequence", f"input:{port}")))
            return None, "", issues
        return model, source, issues

    amber.secondary_structure = secondary_structure


# ---------------------------------------------------------------------------------------------
# script helpers

def _insert_mdbridge(text: str) -> str:
    """MDBRIDGE_LINES after the script header: after `cd "$(dirname "$0")"`, its blank line and the
    environment block (D3, D6)."""
    lines = text.split("\n")
    k = lines.index('cd "$(dirname "$0")"') + 2
    if k < len(lines) and lines[k] == "# Environment":
        while lines[k] != "":
            k += 1
        k += 1
    return "\n".join(lines[:k] + MDBRIDGE_LINES + lines[k:])


def _replace_file(files, name, content):
    import dataclasses

    return [dataclasses.replace(f, content=content) if f.name == name else f for f in files]


def uses_mdbridge(files) -> bool:
    return any(f.role in ("run", "run-separate") and MDBRIDGE_LINES[1] in f.content for f in files)


# ---------------------------------------------------------------------------------------------
# amber-build: WG-06, D3, D9

# MDATK 0.0.5: Mg2+ pair references of the L-RNA library types (Dudek & Trylska style): bridging
# phosphate oxygens O3'/O5', sugar ring O4', hydroxyl oxygens O2' / terminal OH
LRNA_DEFAULT_REFS = {"LD": "OS", "LE": "OS", "LH": "OH"}
CP_ALIASES = {"OP": "O2", "OR": "OS"}
ION_TYPES = ("Mg2+", "MG")


def count_solute_residues(pdb_path):
    """WG-06 (MDATK 0.1.0, without its later virtual-residue option)."""
    from mdatk.workflows.scripts import prep_pdb

    skip = {"HOH", "WAT", "SOL", "NA", "NA+", "CL", "CL-", "K", "K+", "MG", "MG2", "ZN", "CA"} | prep_pdb.WATERS
    count, last = 0, None
    try:
        with open(pdb_path, "r", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if line.startswith("ENDMDL"):
                    break
                if line.startswith("TER"):
                    last = None
                elif line.startswith(("ATOM  ", "HETATM")) and len(line) >= 26:
                    name = line[17:21].strip().upper()
                    if name in skip:
                        continue
                    key = (name, line[21], line[22:27])
                    if key != last:
                        count += 1
                        last = key
    except OSError:
        return None
    return count or None


def _frcmod_text(model, water, water_label, chloride=False):
    """D9 (MDATK 0.0.5 mg_params.frcmod_text): the ion's terms under both type names."""
    from mdatk.workflows import mg_params as M

    p = M.get_set(model, water)
    cl_rmin, cl_eps = p.chloride_lj[0] * 10.0 * M.TWO_SIXTH / 2.0, p.chloride_lj[1] / M.KJ_PER_KCAL
    return (f"{M.MODELS[model]} Mg2+ for {water_label} water ({M.CITATIONS[water]}); "
            f"sigma {p.ion[0]} nm, epsilon {p.ion[1]} kJ/mol"
            + (f"; Cl- {p.chloride}: sigma {p.chloride_lj[0]} nm, epsilon {p.chloride_lj[1]} kJ/mol"
               if chloride else "") + "\n"
            "MASS\n"
            + "".join(f"{t:<4} 24.305\n" for t in ION_TYPES)
            + ("Cl-  35.45\n" if chloride else "") +
            "\n"
            "NONB\n"
            + "".join(f"  {t:<8} {p.rmin_half:.8f}  {p.eps_kcal:.8f}\n" for t in ION_TYPES)
            + (f"  Cl-      {cl_rmin:.8f}  {cl_eps:.8f}\n" if chloride else "") +
            "\n")


def _parmed_mg_lines(original, model, water, water_label, frcmod, chloride, extra_types):
    """D9 (MDATK 0.0.5 mg_params.parmed_lines): the pair terms of the twins of non-standard types."""
    from mdatk.workflows import mg_params as M

    lines = original(model, water, water_label, frcmod, False)
    if extra_types:
        p = M.get_set(model, water)
        lines.append("# non-standard atom types (e.g. an L-RNA library) with the terms of their chemical twin")
        for ty, ref in extra_types.items():
            rmin, e = M.amber_pair(*p.rna[CP_ALIASES.get(ref, ref)])
            lines += [f"# {ty} as {ref}", f"changeLJPair :MG @%{ty.replace('*', chr(92) + '*')} {rmin:.8f} {e:.8f}"]
    if chloride:
        p = M.get_set(model, water)
        rmin, e = M.amber_pair(*p.mg_cl)
        lines += [f"# Mg2+-Cl- pair term (optimised with the {p.chloride} Cl- model)",
                  f"changeLJPair :MG :Cl- {rmin:.8f} {e:.8f}"]
    return lines


def _lrna_tokens(s):
    """D9 (MDATK 0.0.5): (type, reference or '') of the L-RNA types, [] when the L-RNA fixes are off;
    independent of the L-nucleotide naming."""
    if not s["lrna_fixes"]:
        return []
    return [tuple(t.split(":", 1)) if ":" in t else (t, "")
            for t in re.split(r"[\s,]+", s["lrna_types"].strip()) if t]


def _mg_type_map(s):
    """Non-standard atom type -> reference type for the Mg2+ pair terms, and the types without one."""
    out, missing = {}, []
    for t, ref in _lrna_tokens(s):
        ref = ref or LRNA_DEFAULT_REFS.get(t, "")
        if ref:
            out[t] = ref
        else:
            missing.append(t)
    return out, missing


# MDATK 0.0.11 scripts/lrna_names.py (the parsing the planner uses)
_LOAD = re.compile(r'^\s*(loadoff|loadamberprep|source)\s+("[^"]*"|\S+)', re.IGNORECASE | re.MULTILINE)
_ASSIGN = re.compile(r'^\s*([^\s=#{}"]+)\s*=', re.MULTILINE)
_MAP_ENTRY = re.compile(r'\{\s*(?:\d+\s+)?"?([^"\s{}]+)"?\s+"?([^"\s{}]+)"?\s*\}')
LANTHANIDE_IONS = {"LA": "lanthanum", "LU": "lutetium"}


def _leap_code(text):
    return "\n".join(line.split("#", 1)[0] for line in text.splitlines())


def _before_pdb(text, pdb=""):
    code = _leap_code(text)
    found = list(re.finditer(r"^.*\bloadpdb\b.*$", code, re.IGNORECASE | re.MULTILINE))
    named = [m for m in found if pdb and os.path.basename(pdb) in m.group(0)]
    m = (named or found or [None])[0]
    return code[:m.start()] if m else code


def _leap_sources(text, pdb=""):
    code = _before_pdb(text, pdb)
    libs, preps, sourced = [], [], []
    for cmd, arg in _LOAD.findall(code):
        {"loadoff": libs, "loadamberprep": preps, "source": sourced}[cmd.lower()].append(arg.strip('"'))
    return libs, preps, sourced, set(_ASSIGN.findall(code))


def _lib_units_text(text):
    return set(re.findall(r"^!entry\.([^.\s]+)\.unit\.atoms table", text, re.MULTILINE))


def _pdb_res_maps(text):
    out = {}
    code = _leap_code(text)
    for m in re.finditer(r"\baddpdbresmap\b", code, re.IGNORECASE):
        start = code.find("{", m.end())
        if start < 0:
            continue
        depth, end = 0, len(code)
        for k in range(start, len(code)):
            depth += {"{": 1, "}": -1}.get(code[k], 0)
            if depth == 0:
                end = k
                break
        for name, unit in _MAP_ENTRY.findall(code[start + 1:end]):
            out.setdefault(name, set()).add(unit)
    return out


def _unmapped(names, units, maps):
    return [m for m in names if m not in units and not (maps.get(m, set()) & units)]


def _consequence(missing):
    silent = [m for m in missing if m in LANTHANIDE_IONS]
    rest = [m for m in missing if m not in LANTHANIDE_IONS]
    text = "tleap would "
    if silent:
        ions = " and ".join(LANTHANIDE_IONS[m] for m in silent)
        text += (f"silently read {' and '.join(silent)} as the {ions} ion of AmberTools' atomic_ions.lib"
                 + (" and " if rest else ""))
    if rest:
        text += "invent parameterless units for " + ("the rest" if silent else ", ".join(rest))
    return text


def _map_suggestion(missing, units):
    entries = [f'{{ "{m}" "{m[1]}0{m[2:]}" }}' for m in missing if f"{m[1]}0{m[2:]}" in units]
    if entries:
        return f"addPdbResMap {{ {' '.join(entries)} }}"
    return 'addPdbResMap { { "LG" "<the library\'s L-G unit>" } ... }'


def _unit_list(units, limit=12):
    return ", ".join(sorted(units)[:limit]) + (" ..." if len(units) > limit else "")


def _lib_paths(paths, ctx):
    found, unread = [], []
    for path in paths:
        path = path.strip()
        here = os.path.join(ctx.root, ctx.folder, path) if getattr(ctx, "root", "") else ""
        for cand in (path, os.path.join(ctx.folder, path) if getattr(ctx, "folder", "") else "", here):
            if cand and os.path.isfile(cand):
                found.append(cand)
                break
        else:
            unread.append(path)
    return found, unread


def _lib_units(paths, ctx):
    found, unread = _lib_paths(paths, ctx)
    units = set()
    for cand in found:
        try:
            with open(cand, encoding="utf-8", errors="replace") as fh:
                units |= _lib_units_text(fh.read())
        except OSError:
            unread.append(cand)
    return units, unread


def _own_libs(s):
    from mdatk.workflows.base import extra_lines

    libs = _leap_sources(s["leap_before_params"] + "\n" + s["leap_before_load"])[0]
    return [f.strip() for f in extra_lines(s["extra_libs"])] + libs


def _check_lrna_names(s, ctx, path):
    """D9 (MDATK 0.0.11 AmberBuildTool._check_lrna_names, the case of an existing input PDB)."""
    from mdatk import pdbtools
    from mdatk.workflows.schema import Issue

    try:
        named = pdbtools.clean(path, scheme="amber", l_ranges=s["lrna_ranges"] or None,
                               terminal_names=s["lrna_terminals"]).structure
    except (OSError, pdbtools.PDBError):
        return []
    names = sorted({r.name for r in named.residues() if r.chirality == "L"})
    leap = s["leap_before_params"] + "\n" + s["leap_before_load"]
    units, unread = _lib_units(_own_libs(s), ctx)
    units |= _leap_sources(leap)[3]
    missing = _unmapped(names, units, _pdb_res_maps(leap))
    if not missing:
        return []
    if unread:
        return [Issue("warning", f"cannot check that the L residue names ({', '.join(names)}) are units of "
                                 f"{', '.join(unread)} until the file can be read (a name in leap's search path "
                                 "is not read)", ("extra_libs",))]
    fix = f"Add under 'Leap commands before loadpdb': {_map_suggestion(missing, units)}"
    return [Issue("error", f"the L-RNA library does not define {', '.join(missing)}, the names MDATK writes for "
                           f"the L residues, and no addPdbResMap maps them onto its units "
                           f"({_unit_list(units)}). {_consequence(missing)}. {fix}",
                  ("leap_before_load", "extra_libs"))]


def _patch_amber_build():
    from mdatk import pdbtools
    from mdatk.workflows import amber, mg_params
    from mdatk.workflows.base import extra_lines
    from mdatk.workflows.schema import Issue

    A = amber.AmberBuildTool
    amber.count_solute_residues = count_solute_residues  # WG-06

    # D9: split_types on by default
    for section in A.sections:
        for opt in section.options:
            if opt.key == "split_types":
                opt.default = True

    # D9: frcmod with both ion type names, pair terms of the twins
    mg_params.frcmod_text = _frcmod_text
    original_parmed_lines = mg_params.parmed_lines

    A._lrna_tokens = staticmethod(_lrna_tokens)
    A._lrna_types = staticmethod(lambda s: [t for t, _ref in _lrna_tokens(s)])
    A._mg_type_map = staticmethod(_mg_type_map)

    def _parmed(self, s, title, out_parm, wet, ctx=None):
        """MDATK 0.0.4.1's, with 0.0.5.1's dry pass without Mg2+ lines and 0.0.5's twin pair terms."""
        lines = ["#", f"# {title}", f"# ParmEd input written by {amber.PROJECT_ABBREV} {amber.__version__}"]
        cp = self._cp_lib(s)
        if s["nbfix19"] and s["split_types"]:
            lines += ["#", "# Give HO and O2 their own Lennard-Jones types so the pair changes below do not also",
                      "# apply to atom types tleap merged with them (water HW, carbonyl O)",
                      "addLJType @%HO"] + ([] if cp else ["addLJType @%O2"])
        if s["nbfix19"]:
            lines += amber._cp_masks(amber.NBFIX19_LINES) if cp else amber.NBFIX19_LINES
        if s["bph0"]:
            lines += amber._cp_masks(amber.BPH0_LINES) if cp else amber.BPH0_LINES
        if (s["nbfix19"] or s["bph0"]) and self._lrna_types(s):
            lines += amber.lrna_pair_lines(s["nbfix19"], s["bph0"], self._lrna_types(s)[0])
        elif s["nbfix19"] or s["bph0"]:
            lines += amber.LRNA_NOTE
        if wet and self._solvated(s) and s["mg_model"] == "1264":
            lines += ["#", "# 12-6-4 Lennard-Jones terms for Mg2+ (Li & Merz)",
                      f"add12_6_4 :MG watermodel {amber.WATERS[s['water']][3]}"
                      + (" polfile lj_1264_pol_cp.dat" if self._cp_lib(s) else "")]
        if wet and s["hmr"]:
            lines += ["#", "# Hydrogen mass repartitioning (4 fs time steps with SHAKE)", "HMassRepartition"]
        extra = extra_lines(s["parmed_extra"])
        if not wet:  # the dry topology has no ions: Mg2+ lines would select nothing
            extra = [ln for ln in extra if not re.search(r":MG\b", ln, re.IGNORECASE)]
        if extra:
            lines += ["#", "# Extra commands"] + extra
        if wet and self._mg_pairs(s):
            lines += _parmed_mg_lines(original_parmed_lines, s["mg_model"], s["water"],
                                      amber.WATERS[s["water"]][0], self._mg_frcmod(ctx), s["mg_cl_pair"],
                                      self._mg_type_map(s)[0])
        lines += ["#", f"outparm {out_parm}"]
        return "\n".join(lines) + "\n"

    A._parmed = _parmed

    def _check_mg_pairs(self, s):
        """MDATK 0.0.4.1's, with 0.0.5's references of the L-RNA types (D9)."""
        out = []
        model, water = s["mg_model"], s["water"]
        name, wlabel = mg_params.MODELS[model], amber.WATERS[water][0]
        if (model, water) in mg_params.UNAVAILABLE:
            out.append(Issue("error", f"{name} is not available for {wlabel}: "
                                      f"{mg_params.UNAVAILABLE[(model, water)]}", ("mg_model", "water")))
            return out
        if (model, water) not in mg_params.MG_SETS:
            have = sorted({amber.WATERS[w][0] for m, w in mg_params.MG_SETS if m == model})
            out.append(Issue("error", f"no published {name} parameters for {wlabel} (available: "
                                      f"{', '.join(have)})", ("mg_model", "water")))
            return out
        if s["mg_count"] == 0:
            out.append(Issue("warning", f"{name} is selected but no Mg²⁺ is added; it only matters when "
                                        "the input already contains MG ions", ("mg_model", "mg_count")))
        if s["rna_ff"] not in ("ol3", "ol3cp"):
            out.append(Issue("warning", f"the {name} Mg²⁺-RNA pair terms were fitted with parmBSC0χOL3 "
                                        "atom types; RNA atom types of this force field that are not "
                                        "in that set keep the plain combination rule", ("mg_model", "rna_ff")))
        mg_extra = [ln for ln in extra_lines(s["parmed_extra"]) if re.search(r":MG\b|Mg2\+", ln, re.IGNORECASE)]
        if mg_extra:
            out.append(Issue("warning", f"the {name} pair terms are applied after the extra ParmEd commands "
                                        "and replace any Mg²⁺ pair terms set there", ("parmed_extra",)))
        p = mg_params.get_set(model, water)
        has_cl = bool(s["mg_balance"] and s["mg_count"] or s["salt_pairs"])
        if s["mg_cl_pair"] and has_cl:
            out.append(Issue("info", f"{name} for {wlabel} was optimised with the {p.chloride} Cl⁻ model; "
                                     "the Mg²⁺-Cl⁻ pair term is applied, other Cl⁻ interactions use Amber's "
                                     "Cl⁻", ("mg_cl_pair",)))
        if not s["mg_cl_pair"] and has_cl:
            out.append(Issue("warning", f"without the pair term, Mg²⁺-Cl⁻ uses the combination rule with the "
                                        f"optimised {name} parameters - not a published combination (the "
                                        f"authors' pair: σ {p.mg_cl[0]:.3f} nm, ε {p.mg_cl[1]:.2f} kJ/mol) and "
                                        "far too strong Mg²⁺-Cl⁻ contacts", ("mg_cl_pair",)))
        refs = sorted(set(p.rna) | set(CP_ALIASES))
        mapping, missing = _mg_type_map(s)
        bad = [f"{t}:{r}" for t, r in mapping.items() if r not in refs or not re.fullmatch(r"[A-Za-z0-9*+-]{1,4}", t)]
        if bad:
            out.append(Issue("error", f"Mg²⁺ pair types {', '.join(bad)}: write TYPE:REFERENCE with a reference "
                                      f"among {', '.join(refs)}", ("lrna_types",)))
        if missing:
            out.append(Issue("error", f"L-RNA atom type{'s' if len(missing) > 1 else ''} {', '.join(missing)} "
                                      f"{'have' if len(missing) > 1 else 'has'} no Mg²⁺ pair reference: without one "
                                      f"{name} binds those atoms several times too strongly; write e.g. "
                                      f"{missing[0]}:OS under L-RNA atom types", ("lrna_types",)))
        if mapping and not bad:
            out.append(Issue("info", f"Mg²⁺ pair terms for non-standard atom types: "
                                     f"{', '.join(f'{t} as {r}' for t, r in mapping.items())}", ("lrna_types",)))
        if self._cp_lib(s):
            out.append(Issue("info", "phosphate oxygens retyped by the CP library (OP, OR) get the authors' O2 / "
                                     "OS terms - a choice by MDATK, since the authors fitted Mg²⁺ against plain "
                                     "OL3 phosphate oxygens", ("mg_model", "cp_source")))
        if extra_lines(s["extra_libs"]) or extra_lines(s["extra_frcmod"]):
            out.append(Issue("warning", f"the {name} pair terms cover the 27 parmBSC0χOL3 RNA atom types and the "
                                        "L-RNA atom types; other atom types of extra residue libraries or parameter "
                                        "files keep the plain combination rule with Mg²⁺",
                             ("mg_model", "extra_libs")))
        out.append(Issue("info", f"{name} ({mg_params.CITATIONS[water]}): Mg²⁺-water uses the "
                                 "Lorentz-Berthelot rule, Mg²⁺-RNA the authors' pair terms; GROMACS "
                                 "conversion keeps them as [ nonbond_params ]", ("mg_model",)))
        return out

    A._check_mg_pairs = _check_mg_pairs

    def _check_lrna(self, s, ctx):
        """MDATK 0.0.4.1's, with 0.0.5's TYPE:REFERENCE tokens and 0.0.11's lanthanide trap (D9)."""
        out = []
        if not s["lrna"]:
            return out
        try:
            pdbtools.parse_ranges(s["lrna_ranges"] or None)
        except pdbtools.PDBError as exc:
            out.append(Issue("error", f"L ranges: {exc}", ("lrna_ranges",)))
            return out
        libs = extra_lines(s["extra_libs"]) or extra_lines(s["extra_frcmod"])
        path = ctx.existing("structure")
        n_l = 0
        if path:
            try:
                struct = pdbtools.read_pdb(path)
                notes = pdbtools.assign_chirality(struct, detect=True, l_ranges=s["lrna_ranges"] or None)
            except (OSError, pdbtools.PDBError) as exc:
                out.append(Issue("warning", f"could not check the input PDB: {exc}", ("lrna",)))
                return out
            n_l = sum(1 for r in struct.residues() if r.chirality == "L")
            for note in notes:
                if note.startswith("WARNING:"):
                    out.append(Issue("warning", note[len("WARNING:"):].strip(), ("lrna_ranges",)))
            if n_l:
                out.append(Issue("info", f"{os.path.basename(path)}: {notes[0]}", ("lrna",)))
        elif s["lrna_ranges"].strip():
            n_l = 1
        if n_l and not libs:
            out.append(Issue("error" if path else "warning",
                             "L-RNA residues need their library and parameters: add them under Force field "
                             "(e.g. Dudek & Trylska's nucleic12_LRNA.lib with parm10_LRNA.dat, or your own "
                             "LRNA.lib + frcmod.LRNA)", ("extra_libs", "extra_frcmod")))
        elif n_l and path and _own_libs(s):
            out += _check_lrna_names(s, ctx, path)
        if n_l and "LD" not in s["leap_before_params"] and not s["lrna_fixes"]:
            out.append(Issue("info", "L-RNA libraries use the atom types LD / LH; if leap reports unknown types, "
                                     "add  addAtomTypes { { \"LD\" \"O\" \"sp3\" } { \"LH\" \"O\" \"sp3\" } }  "
                                     "under 'Leap commands before parameters'", ("leap_before_params",)))
        bad = [f"{t}:{r}" if r else t for t, r in _lrna_tokens(s)
               if not re.fullmatch(r"[A-Za-z0-9*+-]{1,4}", t) or (r and not re.fullmatch(r"[A-Za-z0-9*+-]{1,4}", r))]
        if bad:
            out.append(Issue("error", f"L-RNA atom types {', '.join(bad)}: 1-4 letters / digits each (TYPE or "
                                      "TYPE:REFERENCE)", ("lrna_types",)))
        elif s["lrna"] and s["lrna_fixes"] and not self._lrna_types(s) and (s["nbfix19"] or s["bph0"]):
            out.append(Issue("warning", "no L-RNA atom type given: the NBfix19 / 0BPh pairs miss L residues",
                             ("lrna_types",)))
        if not path and ctx.has("structure"):
            out.append(Issue("info", "L-nucleotides are detected from the sugar geometry when the script runs "
                                     "(the input PDB does not exist yet)", ("lrna",)))
        return out

    A._check_lrna = _check_lrna

    split_info = ("tleap gives HO the same Lennard-Jones type as water HW (both zero) and O2 the same as carbonyl O, "
                  "so 'changeLJPair @%HO ...' also changes water-H and carbonyl-O pairs; tick 'Separate HO and O2 "
                  "Lennard-Jones types first' if that is not intended")
    split_warning = ("tleap gives HO the same Lennard-Jones type as water HW (both zero) and O2 the same as carbonyl "
                     "O, so the NBfix19 'changeLJPair @%HO ...' lines also apply to every water hydrogen (and "
                     "carbonyl O) for the whole run; tick 'Separate HO and O2 Lennard-Jones types first' unless you "
                     "reproduce a protocol that did not separate them")
    check = A.check

    def check_d9(self, s, ctx):
        out = check(self, s, ctx)
        return [Issue("warning", split_warning, i.keys) if i.message == split_info else i for i in out]

    A.check = check_d9

    generate = A.generate

    def generate_d3(self, s, ctx):
        files = generate(self, s, ctx)
        if not s["lrna"]:
            return files
        n = ctx.name
        tools_py = n("pdbtools", "py")
        script = n(self.tag, "sh")
        files = [f for f in files if f.name != tools_py]
        text = next(f.content for f in files if f.name == script)
        text = _sub(text, f"python3 {tools_py} clean ", f"{MDBRIDGE_WORD} pdb clean ")
        text = _sub(text, "ERROR: pdbtools did not write", "ERROR: mdbridge pdb clean did not write")
        return _replace_file(files, script, _insert_mdbridge(text))

    A.generate = generate_d3


# ---------------------------------------------------------------------------------------------
# amber-md: WG-01 (DISANG), WG-04, WG-05

def _patch_amber_md():
    from mdatk.workflows import amber
    from mdatk.workflows.base import require_file, script_header
    from mdatk.workflows.schema import GeneratedFile, Issue

    T = amber.AmberMDTool
    check = T.check

    def check_fixed(self, s, ctx):
        out = check(self, s, ctx)
        stages = s["stages"]
        if not stages:
            return out
        tags = [self.stage_tag(st, k) for k, st in enumerate(stages)]
        modes = self._restart_modes(s, ctx)
        # WG-05 (MDATK 0.1.0)
        for k, st in enumerate(stages):
            nseg = st["segments"] if st["kind"] == "md" else 1
            if nseg < 2 or tags.count(tags[k]) > 1:
                continue
            parts = [(f"{tags[k]}-{seg:03d}", f"segment {seg}", False) for seg in range(1, nseg + 1)]
            if modes[k][0] == 0 or st["ramp"]:
                parts.append((f"{tags[k]}-cont", f"its input for segments 2-{nseg}", True))
            for tag, what, is_input in parts:
                for j, other in enumerate(stages):
                    one_run = other["kind"] != "md" or other["segments"] < 2
                    if tags[j] == tag and (is_input or one_run):
                        out.append(Issue("error", f"stages {tags[k]} ({what}) and {tag} would write the same files "
                                                  f"({ctx.project}_{ctx.step:02d}_{tag}_…): rename one of them",
                                         (f"stages.{k}.name", f"stages.{j}.name")))
        # WG-01: DISANG= cannot hold white space
        disang = ctx.inp("restraints")
        needs = [k for k, st in enumerate(stages) if st["distres"]]
        if needs and disang and re.search(r"\s", disang):
            out.append(Issue("error", "Restraints (DISANG): pmemd reads DISANG= up to the first space, and this path "
                                      f"has one: {disang}; give the file a path without spaces",
                             tuple(f"stages.{k}.distres" for k in needs) + ("input:restraints",)))
        return out

    T.check = check_fixed

    def generate(self, s, ctx):
        """MDATK 0.0.4.1's AmberMDTool.generate with 0.1.0's continuation input after a ramp (WG-04)."""
        files = []
        runs = self._segments(s, ctx)
        modes = self._restart_modes(s, ctx)
        disang = ctx.inp("restraints")
        solute = self._solute_mask(ctx)
        top, crd = ctx.inp("topology"), ctx.inp("coordinates")
        mdin_of = {}
        for k, st in enumerate(s["stages"]):
            irest, ntx = modes[k]
            tag = self.stage_tag(st, k)
            name = ctx.name(tag, "in")
            files.append(GeneratedFile(name, self.mdin(st, irest, ntx, disang, solute), "input",
                                       f"pmemd input of stage {tag}"))
            nseg = st["segments"] if st["kind"] == "md" else 1
            for seg in range(1, nseg + 1):
                mdin_of[(k, seg)] = name
            if nseg > 1 and (irest == 0 or st["ramp"]):
                cont = ctx.name(f"{tag}-cont", "in")
                files.append(GeneratedFile(cont, self.mdin(dict(st, ramp=False), 1, 5, disang, solute), "input",
                                           f"pmemd input of stage {tag}, segments 2-{nseg}"))
                for seg in range(2, nseg + 1):
                    mdin_of[(k, seg)] = cont
        engine = s["engine"]
        if engine == "pmemd.MPI":
            launcher = s["mpi_launcher"].replace("{n}", str(s["nproc"])).strip()
            pmemd = f"{launcher} pmemd.MPI".strip()
        else:
            pmemd = engine
        min_pmemd = "pmemd" if s["min_engine"] == "pmemd" else "$PMEMD"
        prod_split = s["separate_production"] and any(r["stage"]["production"] for r in runs)

        def preamble(title, script):
            lines = script_header(title, ctx, s, script, [f"Amber {s['amber_version']}; engine {engine}."])
            if engine == "pmemd.cuda" and s["gpu"].strip():
                lines += [f"export CUDA_VISIBLE_DEVICES={s['gpu'].strip()}", ""]
            lines += [f'PMEMD="{pmemd}"', f'TOP="{top}"', ""]
            if s["resume"]:
                lines += ["# a stage counts as finished when its output has the final timing lines",
                          '# (pmemd: "Total wall time", sander: "Run   done at")',
                          'finished() { [ -s "$1" ] && grep -qE "Total wall time|Run +done at" "$1"; }', ""]
            return lines

        main_name = ctx.name(self.tag, "sh")
        prod_name = ctx.name("amber-prod", "sh")
        main = preamble("Amber MD stages", main_name)
        prod = preamble("Amber production", prod_name) if prod_split else main
        prev_rst = crd
        for r in runs:
            st = r["stage"]
            target = prod if (prod_split and st["production"]) else main
            ref = crd if s["restraint_ref"] == "initial" else prev_rst
            exe = min_pmemd if st["kind"] == "min" else "$PMEMD"
            seg = f" (segment {r['seg']}/{r['nseg']})" if r["nseg"] > 1 else ""
            cmd = (f'{exe} -O -i {mdin_of[(r["index"], r["seg"])]} -p "$TOP" -c {prev_rst} -ref {ref} '
                   f'-o {r["out"]} -r {r["rst"]}' + (f" -x {r['nc']}" if r["nc"] else "") + f' -inf {r["inf"]}')
            target.append(f"# {r['rtag']}{seg}: {st['title'] or st['name']}")
            if s["resume"]:
                target += [f"if finished {r['out']}; then", f'  echo "{r["rtag"]}: already finished - skipping"',
                           "else", f"  {cmd}", "fi"]
            else:
                target.append(cmd)
            target += [require_file(r["rst"], f"pmemd ({r['rtag']})"), ""]
            prev_rst = r["rst"]
        main_runs = [r for r in runs if not (prod_split and r["stage"]["production"])]
        if main_runs and s["dry_pdb"]:
            dry = ctx.name("dry", "pdb")
            main += [f"# Strip a dry PDB of {main_runs[-1]['rtag']} for visual checks",
                     "cpptraj << EOF", 'parm "$TOP"', f"trajin {main_runs[-1]['rst']}",
                     f"strip {re.sub(r'[ ]+', '', s['strip_mask'])}", f"trajout {dry} pdb", "run",
                     "EOF", require_file(dry, "cpptraj"), ""]
        prod_runs = [r for r in runs if r["nc"] and r["stage"]["production"]]
        if s["dry_traj"] and prod_runs:
            last = prod if prod_split else main
            mask_ = re.sub(r"[ ]+", "", s["dry_mask"])
            last += ["# Dry, imaged copy of the production trajectory (and its topology)",
                     "cpptraj << EOF", 'parm "$TOP"'] + [f"trajin {r['nc']}" for r in prod_runs] + [
                "autoimage", f"strip {mask_} parmout {ctx.name('dry', 'parm7')}",
                f"trajout {ctx.name('dry-prod', 'nc')} netcdf", "run", "EOF",
                require_file(ctx.name("dry-prod", "nc"), "cpptraj"), ""]
        files.append(GeneratedFile(main_name, "\n".join(main), "run", "runs the MD stages"))
        if prod_split:
            files.append(GeneratedFile(prod_name, "\n".join(prod), "run-separate", "runs the production stages"))
        return files

    T.generate = generate


# ---------------------------------------------------------------------------------------------
# FARFAR2, FARFAR2 models, SimRNA, RNAMake, VMD: WG-07, WG-08, WG-09

def _failed_runs(program, runs, fewer, missing):
    return [f"# a failed run does not stop the others and leaves fewer {fewer}; the step fails when every run failed",
            None,  # the xargs line
            'NFAILED=$(ls "$FAILED" | wc -l); rm -rf "$FAILED"',
            f'if [ "$NFAILED" -ge {runs} ]; then echo "ERROR: every {program} run failed (see above)" >&2; exit 1; fi',
            f'if [ "$NFAILED" -gt 0 ]; then echo "WARNING: $NFAILED of {runs} {program} runs failed (see above): '
            f'{missing} are missing" >&2; fi']


def _xargs_failures(text, s, program, fewer, missing, exports_old, exports_new):
    xargs = f"seq 1 {s['runs']} | xargs -P {s['parallel']} -I{{}} bash -c 'run_one {{}}'"
    block = _failed_runs(program, s["runs"], fewer, missing)
    block[1] = xargs + " || true"
    text = _sub(text, f"\n{exports_old}\n", f"\nFAILED=$(mktemp -d)\n{exports_new}\n")
    return _sub(text, f"\n{xargs}\n", "\n" + "\n".join(block) + "\n")


POSRE_OLD = """    if len(core) > 1 and core[0] in "RDL":
        core = core[1:]
"""
POSRE_NEW = """    if len(core) > 2 and core.startswith("LD"):  # L-DNA: LDA LDC LDG LDT (the units of MDATK's L-DNA library)
        core = core[2:]
    elif len(core) > 1 and core[0] in "RDL":
        core = core[1:]
"""


def _patch_modelling():
    from mdatk.workflows import farfar2, gromacs, rnamake, simrna, vmd

    G = gromacs.GromacsTool if hasattr(gromacs, "GromacsTool") else next(
        v for v in vars(gromacs).values() if isinstance(v, type) and getattr(v, "name", "") == "gromacs")
    generate_g = G.generate

    def generate_gromacs(self, s, ctx):  # WG-21 (MDATK 0.1.0)
        files = generate_g(self, s, ctx)
        helper = ctx.name("posre", "py")
        return [_replace_file([f], helper, _sub(f.content, POSRE_OLD, POSRE_NEW))[0] if f.name == helper else f
                for f in files]

    G.generate = generate_gromacs

    F = farfar2.FarfarTool
    generate_f = F.generate

    def generate_farfar2(self, s, ctx):  # WG-07 (MDATK 0.1.0)
        files = generate_f(self, s, ctx)
        script = ctx.name(self.tag, "sh")
        text = next(f.content for f in files if f.name == script)
        text = _sub(text, '|| return 0\n', '|| return 1\n')
        log = re.search(r'\|\| echo "run \$N failed - see ([^"]+)" >&2\n', text)
        text = _sub(text, log.group(0), f'|| {{ echo "run $N failed - see {log.group(1)}" >&2; : > "$FAILED/$N"; '
                                        'return 1; }\n')
        text = _xargs_failures(text, s, "rna_denovo", "models", "their models", "export FLAGS", "export FLAGS FAILED")
        return _replace_file(files, script, text)

    F.generate = generate_farfar2

    M = farfar2.FarfarModelsTool
    generate_m = M.generate

    def generate_models(self, s, ctx):  # WG-09 (MDATK 0.1.0)
        files = generate_m(self, s, ctx)
        script = ctx.name(self.tag, "sh")
        text = next(f.content for f in files if f.name == script)
        every = ctx.name("model-MODELNUM", "pdb").replace("MODELNUM", "[0-9]*")
        old = f"# Extract the best {s['n_models']} models as numbered PDB files\nK=0\n"
        new = (f"# Extract the best {s['n_models']} models as numbered PDB files (those of an earlier run go first)\n"
               f"rm -f {ctx.name('models', 'txt')} {every}\nK=0\n")
        text = _sub(text, old, new)
        if s["dssr"]:  # WG-19 (MDATK 0.1.0)
            pattern = ctx.name("model-MODELNUM", "pdb").replace("MODELNUM", "[0-9][0-9]")
            text = _sub(text, f"# Annotate every model with DSSR\nfor f in {pattern}; do\n",
                        f"# Annotate every model with DSSR\nfor f in {every}; do\n")
        return _replace_file(files, script, text)

    M.generate = generate_models

    S = simrna.SimRNATool
    generate_s = S.generate

    def generate_simrna(self, s, ctx):  # WG-07, WG-08, WG-09 (MDATK 0.1.0)
        files = generate_s(self, s, ctx)
        script = ctx.name(self.tag, "sh")
        text = next(f.content for f in files if f.name == script)
        text = _sub(text, '|| return 0\n', '|| return 1\n')
        log = re.search(r'\|\| echo "run \$N failed - see ([^"]+)" >&2\n', text)
        text = _sub(text, log.group(0), f'|| {{ echo "run $N failed - see {log.group(1)}" >&2; : > "$FAILED/$N"; '
                                        'return 1; }\n')
        text = _xargs_failures(text, s, "SimRNA", "frames", "their trajectories", "export SIMRNA_DIR",
                               "export SIMRNA_DIR FAILED")
        if not s["cluster"]:  # WG-20 (MDATK 0.1.0): the output 'trajectory' (with clustering its script writes it)
            from mdatk.workflows.base import require_file

            all_trafl = ctx.name("simrna-all", "trafl")
            text += "\n".join(["", "# all replica trajectories of all runs", f"cat runs/r*/*.trafl > {all_trafl}",
                               require_file(all_trafl, "SimRNA"), ""])
        files = _replace_file(files, script, text)
        if s["cluster"]:
            cluster = ctx.name("simrna-cluster", "sh")
            text = next(f.content for f in files if f.name == cluster)
            text = _sub(text, "REF=$(ls runs/r*/*-000001.pdb 2>/dev/null | head -n 1)\n",
                        "REF=$(ls runs/r*/*-000001.pdb 2>/dev/null | head -n 1 || true)\n")
            all_trafl = ctx.name("simrna-all", "trafl")
            models = ctx.name("model-MODELNUM", "pdb").replace("MODELNUM", "[0-9]*")
            text = _sub(text, '[ -n "$REF" ] || { echo "ERROR: no first-frame PDB found in runs/" >&2; exit 1; }\n\n',
                        '[ -n "$REF" ] || { echo "ERROR: no first-frame PDB found in runs/" >&2; exit 1; }\n\n'
                        "# the clusters and models of an earlier run go first: another threshold names its clusters\n"
                        "# differently, and the loop below would read them\n"
                        f"rm -f {all_trafl[:-6]}_thrs*_clust* {models}\n")
            files = _replace_file(files, cluster, text)
        return files

    S.generate = generate_simrna

    R = rnamake.RNAMakeTool
    generate_r = R.generate

    def generate_rnamake(self, s, ctx):  # WG-09 (MDATK 0.1.0)
        files = generate_r(self, s, ctx)
        if not s["dump_pdbs"]:
            return files
        script = ctx.name(self.tag, "sh")
        text = next(f.content for f in files if f.name == script)
        pattern = ctx.name("design-MODELNUM", "pdb").replace("MODELNUM", "[0-9]*")
        old = "\n\n# Build motif paths between the base pairs, optimise sequences, score designs\n"
        new = ("\n\n# The design PDBs of an earlier run go first\nrm -f design.*.pdb " + pattern + "\n" + old[1:])
        return _replace_file(files, script, _sub(text, old, new))

    R.generate = generate_rnamake

    V = vmd.VMDTool
    generate_v = V.generate

    def generate_vmd(self, s, ctx):  # WG-09 (MDATK 0.1.0)
        files = generate_v(self, s, ctx)
        if s["render"] != "movie":
            return files
        script = ctx.name(self.tag, "sh")
        text = next(f.content for f in files if f.name == script)
        exe = s["vmd_exe"].strip() or "vmd"
        batch = "-dispdev text " if s["renderer"] == "TachyonInternal" else ""
        line = f"\n{exe} {batch}-e {ctx.name(self.tag, 'tcl')}\n"
        return _replace_file(files, script, _sub(text, line, "\n# the frames of an earlier run go first\n"
                                                             "rm -f frames/frame.*.tga" + line))

    V.generate = generate_vmd


# ---------------------------------------------------------------------------------------------
# structure builder: WG-18 / D6

def _patch_builder():
    from mdatk.workflows import buildstep

    B = buildstep.StructureBuilderTool if hasattr(buildstep, "StructureBuilderTool") else None
    if B is None:
        B = next(v for v in vars(buildstep).values() if isinstance(v, type) and getattr(v, "name", "") ==
                 "structure-builder")

    def runs_mdbridge(s):
        return (s["python"].strip() or "python3") == "python3"

    mdatk_words = B.mdatk

    def mdatk(s):
        return [MDBRIDGE_WORD] if runs_mdbridge(s) else mdatk_words(s)

    B.mdatk = staticmethod(mdatk)
    generate = B.generate

    def generate_d6(self, s, ctx):
        files = generate(self, s, ctx)
        if not runs_mdbridge(s):
            return files
        script = ctx.name(self.tag, "sh")
        text = next(f.content for f in files if f.name == script)
        text = _sub(text, "\n# MDATK structure builder: ", f"\n# {NAME} structure builder: ")
        text = _sub(text, "# MDATK must be importable by the Python below (installed, or its folder on PYTHONPATH).\n",
                    "# MDBridge must be on PATH, or MDBRIDGE must name the program.\n")
        text = _sub(text, "ERROR: mdatk build did not write", "ERROR: mdbridge build did not write")
        return _replace_file(files, script, _insert_mdbridge(text))

    B.generate = generate_d6


# ---------------------------------------------------------------------------------------------
# restraint dictionaries: D3

def _patch_dicts():
    from mdatk.workflows import dictstep

    D = next(v for v in vars(dictstep).values() if isinstance(v, type) and getattr(v, "name", "") == "restraint-dicts")
    generate = D.generate

    def generate_d3(self, s, ctx):
        files = generate(self, s, ctx)
        helper = ctx.name("dictionaries", "py")
        script = ctx.name(self.tag, "sh")
        files = [f for f in files if f.name != helper]
        text = next(f.content for f in files if f.name == script)
        text = _sub(text, f"python3 {helper} finish ", f"{MDBRIDGE_WORD} dict finish ")
        text = _sub(text, "ERROR: dictionaries.py finish did not write", "ERROR: mdbridge dict finish did not write")
        return _replace_file(files, script, _insert_mdbridge(text))

    D.generate = generate_d3


# ---------------------------------------------------------------------------------------------
# the plan: WG-01 (quoting), WG-02 (master), WG-03, D3 / D6 (summary)

_WORD_END = r"(?=[\s;|&)]|$)"


def quote_chosen(text: str, paths) -> str:
    """WG-01: hand-picked paths as words of shell lines, through shell_quote (lines that are not comments,
    outside here-documents; a path already inside quotes or a word of its own keeps its text)."""
    from mdatk.workflows.base import shell_quote

    paths = sorted({p for p in paths if p and shell_quote(p) != p}, key=len, reverse=True)
    if not paths:
        return text
    out, heredoc = [], None
    for line in text.split("\n"):
        if heredoc is not None:
            if line.strip() == heredoc:
                heredoc = None
            out.append(line)
            continue
        m = re.search(r"<<-?\s*'?\"?([A-Za-z_][A-Za-z0-9_]*)", line)
        if not line.lstrip().startswith("#"):
            for p in paths:
                line = re.sub(r"(?<=[\s])" + re.escape(p) + _WORD_END, lambda _m, q=shell_quote(p): q, line)
        if m:
            heredoc = m.group(1)
        out.append(line)
    return "\n".join(out)


def _chosen_paths(step):
    return [p for v in step.ctx.inputs.values() if not v.source for p in v.paths] if step.ctx else []


def broken_sources(plan, nid):
    """The steps with errors a step reads from, directly or through other steps, in plan order (WG-02)."""
    sources = {}
    for e in plan.workflow.edges:
        sources.setdefault(e.target, set()).add(e.source)
    seen, todo = set(), list(sources.get(nid, ()))
    while todo:
        cur = todo.pop()
        if cur in seen or cur not in plan.steps:
            continue
        seen.add(cur)
        todo.extend(sources.get(cur, ()))
    order = plan.order + [n for n in plan.steps if n not in plan.order]
    return [n for n in order if n in seen and plan.steps[n].tool is not None and plan.steps[n].tool.generates
            and not plan.steps[n].ok]


def _patch_plan():
    import dataclasses

    from mdatk.workflows import plan as P
    from mdatk.workflows.schema import Issue

    generate_step = P.generate_step

    def generate_step_quoted(step):  # WG-01
        files = generate_step(step)
        chosen = _chosen_paths(step)
        return [dataclasses.replace(f, content=quote_chosen(f.content, chosen))
                if f.role in ("run", "run-separate") else f for f in files]

    P.generate_step = generate_step_quoted

    check_names = P.Plan._check_steps_and_names

    def _check_steps_and_names(self):  # WG-03 (MDATK 0.1.2 _check_shared_outputs)
        check_names(self)
        claimed = {}
        for nid in self.order:
            step = self.steps[nid]
            if step.tool is None or not step.tool.generates or step.ctx is None:
                continue
            passed = {p for v in step.ctx.inputs.values() for p in v.paths}
            passed |= {str(step.settings.get(key) or "").strip()
                       for key, opt in step.tool.flat_options().items() if opt.kind == "file"}
            for port, value in step.outputs.items():
                for path in value.paths:
                    if os.path.isabs(path) or path in passed or os.path.dirname(path) != step.ctx.folder:
                        continue
                    first = self.steps[claimed.setdefault(path, nid)]
                    if first is step or (first.node.step, first.tool.tag) == (step.node.step, step.tool.tag):
                        continue
                    step.issues.append(Issue("error", f"output '{port}': {path} is also written by "
                                                      f"{first.label}; give one of them another file name"))

    P.Plan._check_steps_and_names = _check_steps_and_names

    def master_script(self):  # WG-02
        lines = ["#!/bin/bash",
                 f"# Master script of plan {self.project} - written by {P.PROJECT_ABBREV} {P.__version__} "
                 f"on {P._dt.date.today().isoformat()}.",
                 "# Runs the main script of every step in plan order.  Scripts marked 'run separately'",
                 "# (e.g. production runs for a cluster) are listed but not started.",
                 "", "set -eo pipefail", 'cd "$(dirname "$0")"', ""]
        for nid in self.order:
            step = self.steps[nid]
            if step.tool is None or not step.tool.generates:
                continue
            lines.append(f"# {step.label}")
            if not step.files and step.ok:
                step.files = P.generate_step(step)
            if not step.ok:
                lines.append("#   (not written: the step has errors)")
            broken = broken_sources(self, nid) if step.ok else []
            if broken:
                labels = [self.steps[n].label for n in broken]
                lines.append(f"#   (not run: reads from {', '.join(labels)}, which "
                             f"{'has' if len(labels) == 1 else 'have'} errors)")
                lines.append("")
                continue
            for f in step.files:
                rel = step.ctx.from_own(f.name)
                if f.role == "run":
                    lines.append(f'bash "{rel}"')
                elif f.role == "run-separate":
                    lines.append(f"#   run separately: bash {rel}")
            lines.append("")
        return "\n".join(lines)

    P.Plan.master_script = master_script

    summary_text = P.Plan.summary_text

    def summary_with_program(self):  # D3 / D6 (WG-18)
        text = summary_text(self)
        lines = text.split("\n")
        for nid in self.order:
            step = self.steps[nid]
            if step.tool is None or not step.files or not uses_mdbridge(step.files):
                continue
            head = lines.index(f"## {step.label}")
            k = lines.index(f"Type: {step.tool.label}", head)
            lines[k:k + 1] = [lines[k] + "  ", MDBRIDGE_SUMMARY]
        return "\n".join(lines)

    P.Plan.summary_text = summary_with_program


# ---------------------------------------------------------------------------------------------
# WG-15: the planner's New stamp

def new_stamp_after(now, current: str) -> str:
    """WG-15: a new stamp never equals the current one: the stamp of `now`, or, when that is the current
    stamp (or earlier), the minute after the current one."""
    from mdatk.workflows.naming import STAMP_FORMAT

    stamp = now.strftime(STAMP_FORMAT)
    if current:
        try:
            cur = _real_dt.datetime.strptime(current, STAMP_FORMAT)
        except ValueError:
            return stamp
        if stamp <= current:
            return (cur + _real_dt.timedelta(minutes=1)).strftime(STAMP_FORMAT)
    return stamp


# ---------------------------------------------------------------------------------------------
# the commands: WG-02 (run), D14 / WG-13

STALE = "its files on disk were written with other settings; write the plan first"
UNWRITTEN = "its files are not written yet; write the plan first"


def stale_steps(plan, ids):
    """D14: (node id, message) of the steps in `ids` (None: all) whose generated files differ from the disk
    (dates and versions masked) or are missing there."""
    out = []
    for nid in plan.order:
        if ids is not None and nid not in ids:
            continue
        step = plan.steps[nid]
        if step.tool is None or not step.tool.generates or not step.ok:
            continue
        if not step.files:
            from mdatk.workflows import plan as P

            step.files = P.generate_step(step)
        why = ""
        for f in step.files:
            path = os.path.join(plan.root, step.ctx.from_own(f.name))
            try:
                with open(path, encoding="utf-8", errors="surrogateescape", newline="") as fh:
                    disk = fh.read()
            except OSError:
                why = why or UNWRITTEN
                continue
            want = f.content if f.content.endswith("\n") else f.content + "\n"
            if mask(disk) != mask(want):
                why = STALE
        if why:
            out.append((nid, f"step {step.node.step:02d}: {why}"))
    return out


def apply_cli():
    global _applied_cli
    apply_mdatk()
    if _applied_cli:
        return
    _applied_cli = True
    from mdatk import cli

    def _run_plan_scripts(plan, ids, project_path, force=False):
        """MDATK 0.0.4.1's, with 0.1.0's handling of steps with errors (WG-02) and D14's stale check."""
        from mdatk.workflows.runner import ScriptFailed, plan_jobs, run_jobs, separate_scripts

        chosen = [n for n in list(plan.order) + [n for n in plan.steps if n not in plan.order]
                  if ids is None or n in ids]
        broken = [n for n in chosen if plan.steps[n].tool is not None and plan.steps[n].tool.generates
                  and not plan.steps[n].ok]
        if broken:
            readers = {}
            for e in plan.workflow.edges:
                readers.setdefault(e.source, set()).add(e.target)
            after, todo = set(), list(broken)
            while todo:
                for nid in readers.get(todo.pop(), ()):
                    if nid not in after:
                        after.add(nid)
                        todo.append(nid)
            for nid in broken:
                why = "; ".join(i.message for i in plan.steps[nid].issues if i.level == "error")
                print(f"error: {plan.steps[nid].label} has errors and was not run: {why}", file=sys.stderr)
            for nid in chosen:
                step = plan.steps[nid]
                if nid in after and nid not in broken and step.tool is not None and step.tool.generates:
                    print(f"{step.label}: not started (it reads from a step with errors)", file=sys.stderr)
            ids = [n for n in chosen if n not in broken and n not in after]
        if not force:
            stale = stale_steps(plan, ids)
            if stale:
                for _nid, message in stale:
                    print(f"error: {message}", file=sys.stderr)
                print("error: nothing was run (use --force to run the scripts on disk anyway)", file=sys.stderr)
                return 1
        jobs = plan_jobs(plan, ids)
        if not jobs:
            if broken:
                return 1
            print("no scripts to run", file=sys.stderr)
            return 1
        try:
            run_jobs(jobs, project_path, on_line=lambda line: print(line, flush=True))
        except ScriptFailed as exc:
            print(f"error: {exc}", file=sys.stderr)
            return 1
        for path in separate_scripts(plan, ids):
            print(f"not started (run separately): {path}")
        return 1 if broken else 0

    cli._run_plan_scripts = _run_plan_scripts

    def cmd_workflow_run(args):
        from mdatk.workflows.plan import Plan

        project = cli._workflow_project(args)
        plan = Plan(project.workflow, project.workflow_dir, project.groups)
        ids = None
        if args.step is not None:
            ids = [n.id for n in project.workflow.nodes if n.step == args.step]
            if not ids:
                raise ValueError(f"no step number {args.step} in the plan")
        return _run_plan_scripts(plan, ids, cli._project_file(args), force=getattr(args, "force", False))

    cli.cmd_workflow_run = cmd_workflow_run
    build_parser = cli.build_parser

    def build_parser_force():
        import argparse

        parser = build_parser()
        for action in parser._actions:
            if isinstance(action, argparse._SubParsersAction) and "workflow" in action.choices:
                for sub in action.choices["workflow"]._actions:
                    if isinstance(sub, argparse._SubParsersAction) and "run" in sub.choices:
                        run = sub.choices["run"]
                        run.add_argument("--force", action="store_true",
                                         help="run the scripts on disk even when the settings would write other "
                                              "files")
                        run.set_defaults(func=cmd_workflow_run)
        return parser

    cli.build_parser = build_parser_force
