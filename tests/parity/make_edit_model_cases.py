"""What MDATK 0.1.1 does where MDBridge's edit model (src/core/project_document.cpp, problems.cpp) follows
it: residue-group renames and references, located problems and their descriptions.

    python3 -I -B tests/parity/make_edit_model_cases.py MDATK_LATEST_CHECKOUT VALIDATE_CASES_JSON OUT_DIR

MDATK_LATEST_CHECKOUT is MDATK at version 0.1.1 (the commit regenerate_all.py names); nothing is written
into it. VALIDATE_CASES_JSON is tests/data/project-model/validate_cases.json (make_model_cases.py), whose
inputs are checked again here. Needs numpy and scipy (0.1.1's density checks import them); MDAnalysis is not
needed:
0.1.1's selection_problems (selection parsing, which needs MDAnalysis) is left out, since MDBridge keeps
v0.0.4.1's blocking rules (owner decision D7 of milestone 5) and compares only the problems both have.
Writes OUT_DIR/edit_model_cases.json:

* ``rename``      Project.rename_group(old, new) on Project.from_dict(input): the changed locations (or
                  the error: KeyError / ValueError and its message) and to_dict() afterwards
* ``references``  Project.group_references(name)
* ``problems``    per validate case: Project.problems() (field, message, location, level) in a folder
                  holding the case's files
* ``describe``    problems.describe(field, project) for fields with and without a project

Tagged values as in make_model_cases.py (tests/core/test_support.hpp).
"""

from __future__ import annotations

import datetime as dt
import json
import os
import sys
import tempfile

if len(sys.argv) != 4:
    sys.exit("usage: python3 -I -B tests/parity/make_edit_model_cases.py MDATK_LATEST VALIDATE_CASES_JSON OUT_DIR")
MDATK, CASES, OUT = (os.path.abspath(a) for a in sys.argv[1:4])
sys.dont_write_bytecode = True
sys.path.insert(0, MDATK)

import mdatk  # noqa: E402
from mdatk import config as C  # noqa: E402
from mdatk import problems as P  # noqa: E402

if mdatk.__version__ != "0.1.1":
    sys.exit(f"{MDATK} is MDATK {mdatk.__version__}, not 0.1.1")

C.Project.selection_problems = lambda self, *args, **kwargs: []  # see the module text

os.makedirs(OUT, exist_ok=True)


def enc(v):
    if v is None:
        return {"t": "none"}
    if isinstance(v, bool):
        return {"t": "bool", "v": v}
    if isinstance(v, int):
        return {"t": "int", "v": str(v)}
    if isinstance(v, float):
        return {"t": "float", "hex": v.hex(), "repr": repr(v)}
    if isinstance(v, str):
        return {"t": "str", "utf8": v.encode("utf-8").hex(), "v": v}
    if isinstance(v, dt.datetime):
        return {"t": "datetime", "iso": v.isoformat()}
    if isinstance(v, (list, tuple)):
        return {"t": "array", "v": [enc(x) for x in v]}
    if isinstance(v, dict):
        return {"t": "table", "v": [[k, enc(x)] for k, x in v.items()]}
    raise TypeError(type(v))


def dec(v):
    t = v["t"]
    if t == "none":
        return None
    if t in ("bool",):
        return v["v"]
    if t == "int":
        return int(v["v"])
    if t == "float":
        return float.fromhex(v["hex"])
    if t == "str":
        return bytes.fromhex(v["utf8"]).decode("utf-8", "surrogateescape")
    if t == "datetime":
        return dt.datetime.fromisoformat(v["iso"])
    if t == "date":
        return dt.date.fromisoformat(v["iso"])
    if t == "time":
        return dt.time.fromisoformat(v["iso"])
    if t == "array":
        return [dec(x) for x in v["v"]]
    if t == "table":
        return {k: dec(x) for k, x in v["v"]}
    raise ValueError(t)


def project_a():
    """A project referring to its groups in every way rename_group follows (and some it must not)."""
    p = C.Project()
    p.notes = "group:Aptamer stays in the notes"
    p.groups = [C.ResidueGroup("Aptamer", "1-12", "", "#2a78d6"),
                C.ResidueGroup("Target", "13-24", "group:Aptamer or group:Aptamerx", "#eb6834"),
                C.ResidueGroup("Aptamer2", "", "group:Aptamer and name P", "#123456"),
                C.ResidueGroup("Ion", "", "resname MG", "#008300")]
    a = p.analyses
    a.rmsd.targets[0].fit_selection = "group:Aptamer2 or group:Aptamer"
    a.rmsd.targets.append(C.RMSDTarget("Bare", "Aptamer", ""))
    a.hbonds.group_b = " Aptamer "
    a.mg.ligand_selection = "Aptamer, Target"
    a.mg.occupancy_selection = "Aptamer"
    a.mg.group_a = "Aptamer2"
    a.stacking.selection = "group:Aptamer"
    a.energy.group_b = "AptamerX"
    v = p.visual
    v.structure_groups = ["Aptamer", "Target", "Aptamer"]
    v.representations = [C.RepSpec("cartoon", ["Aptamer", "Ion"]), C.RepSpec("licorice", ["Target"])]
    v.proximity_set = "Aptamer"
    v.label_residues = "group:Aptamer"
    v.center_on = "Aptamer"
    p.layout.panels = [
        C.PanelSpec("structure", [], "Structure"),
        C.PanelSpec("timeseries", ["contacts@Aptamer+Target/atom_pairs", "rmsd/Complex",
                                   "hbonds@Aptamer2+Aptamer/count", "mg@Aptamer"], "Contacts"),
        C.PanelSpec("contact_map", ["contacts/map"], "",
                    {"mode": "instant", "note": "group:Aptamer", "list": ["group:Aptamer", 1, "Aptamer"],
                     "nested": {"sel": "group:Aptamer and name P"}}),
    ]
    p.output.charts.series = ["contacts@Aptamer+Target/per_residue_a", "energy@Target+Aptamer"]
    p.output.pdb.selection = "Aptamer,Target, Ion"
    p.workflow.nodes = [C.PlanNode(id="n1", tool="note", step=1, settings={"text": "group:Aptamer", "n": 3})]
    p.workflow.scratch = [C.PlanNode(id="s1", tool="note", step=1, adopted="group:Aptamer")]
    return p


def project_b():
    """Selections that are no list of group names and keys that name other groups."""
    p = C.Project()
    p.groups.append(C.ResidueGroup("A.b-1", "1", "", "#000000"))
    a = p.analyses
    a.mg.ligand_selection = "Aptamer, resid 5"
    a.mg.occupancy_selection = "Aptamer,,Target,"
    a.stacking.selection = "group:Aptamer'"
    p.output.pdb.selection = "Aptamer and name P"
    p.visual.proximity_set = "sel:group:Aptamer"
    p.layout.panels[1].series = ["rmsd@Aptamer/Complex", "@Aptamer/x", "a/b@Aptamer", "contacts@A.b-1+Aptamer"]
    return p


def rename_case(name, project, old, new):
    data = project.to_dict()
    p = C.Project.from_dict(data)
    case = {"name": name, "input": enc(data), "old": old, "new": new}
    try:
        case["changed"] = p.rename_group(old, new)
        case["output"] = enc(p.to_dict())
    except KeyError as exc:
        case["error"] = {"type": "KeyError", "message": exc.args[0]}
    except ValueError as exc:
        case["error"] = {"type": "ValueError", "message": str(exc)}
    return case


def duplicate_project():
    p = C.Project()
    p.groups.append(C.ResidueGroup("Aptamer", "30", "", "#000000"))
    return p


RENAME = [
    rename_case("everything", project_a(), "Aptamer", "Apt"),
    rename_case("stripped new name", project_a(), "Aptamer", "  Apt_2  "),
    rename_case("target", project_a(), "Target", "T.1"),
    rename_case("group used by another group", project_a(), "Aptamer2", "Second"),
    rename_case("no list of names", project_b(), "Aptamer", "Apt"),
    rename_case("dotted name", project_b(), "A.b-1", "Ab"),
    rename_case("same name", project_a(), "Aptamer", "Aptamer"),
    rename_case("same name with spaces", project_a(), "Aptamer", " Aptamer "),
    rename_case("unknown group", project_a(), "Nope", "X"),
    rename_case("empty new name", project_a(), "Aptamer", "   "),
    rename_case("invalid new name", project_a(), "Aptamer", "a b"),
    rename_case("taken name", project_a(), "Aptamer", "Target"),
    rename_case("several groups of the old name", duplicate_project(), "Aptamer", "X"),
    rename_case("default project", C.Project(), "Target", "Protein"),
]


def references_case(name, project, group):
    data = project.to_dict()
    return {"name": name, "input": enc(data), "group": group,
            "found": C.Project.from_dict(data).group_references(group)}


def deleted_ion():
    p = project_a()
    p.groups = [g for g in p.groups if g.name != "Ion"]
    return p


REFERENCES = [
    references_case("everything", project_a(), "Aptamer"),
    references_case("target", project_a(), "Target"),
    references_case("deleted group", deleted_ion(), "Ion"),
    references_case("no list of names", project_b(), "Aptamer"),
    references_case("unused", project_a(), "Nobody"),
]


def problems_cases():
    with open(CASES, encoding="utf-8") as fh:
        vectors = json.load(fh)
    out = []
    here = os.getcwd()
    with tempfile.TemporaryDirectory() as tmp:
        for name in vectors["files"]:
            open(os.path.join(tmp, name), "w").close()
        for name in vectors["folders"]:
            os.mkdir(os.path.join(tmp, name))
        os.chdir(tmp)
        try:
            for case in vectors["cases"]:
                entry = {"name": case["name"]}
                try:
                    p = C.Project.from_dict(dec(case["input"]))
                    entry["problems"] = [{"field": x.field, "message": x.message, "location": x.location,
                                          "level": x.level} for x in p.problems()]
                except (ValueError, TypeError, KeyError) as exc:  # 0.1.1 rejects an input 0.0.4.1 took
                    entry["error"] = f"{type(exc).__name__}: {exc}"
                out.append(entry)
        finally:
            os.chdir(here)
    return out


DESCRIBE_FIELDS = [
    "", "notes", "inputs.topology", "inputs.trajectories", "inputs.stop", "inputs.dt_ps", "inputs.external.0",
    "inputs.external.0.path", "cleaning", "cleaning.fit_mask", "groups", "groups.0", "groups.1.name",
    "groups.1.residues", "groups.1.selection", "groups.0.color", "groups.7.name", "analyses.rmsd",
    "analyses.rmsd.targets.0.selection", "analyses.rmsd.targets.2.fit_selection", "analyses.rmsd.targets.1.label",
    "analyses.contacts.cutoff", "analyses.contacts.group_a", "analyses.hbonds.group_b",
    "analyses.hbonds.angle_cutoff", "analyses.mg.inner_cutoff", "analyses.mg.occupancy_selection",
    "analyses.mg.mg_selection", "analyses.stacking.selection", "analyses.basepairs.engine",
    "analyses.energy.dielectric", "analyses.proximity_search_radius", "analyses.use_pbc",
    "visual.representations.0.style", "visual.representations.1.groups", "visual.representations.0.color",
    "visual.representations.0.size", "visual.structure_groups", "visual.proximity_set",
    "visual.proximity_cutoff", "visual.color_mode", "visual.zoom", "visual.renderer", "visual.rotation",
    "density", "density.map_path", "density.transform", "density.fit_frame", "layout.preset", "layout.panels",
    "layout.panels.0", "layout.panels.1.series", "layout.panels.2.options", "layout.panels.2.options.mode",
    "layout.panels.2.groups", "layout.panels.9.title", "output.video.fps", "output.video.width",
    "output.charts.formats", "output.pdb.selection", "output.out_dir", "workflow.nodes.0.settings.text",
    "workflow.parallel", "nonsense.deeper_key",
]


def describe_cases():
    p = project_a()
    out = []
    for field in DESCRIBE_FIELDS:
        out.append({"field": field, "plain": P.describe(field), "project": P.describe(field, p)})
    return out


def dump(name, data):
    lines = []
    for key, value in data.items():
        head = json.dumps(key, ensure_ascii=False) + ": "
        if isinstance(value, list) and value:
            items = [" " + json.dumps(v, ensure_ascii=False) for v in value]
            lines.append(head + "[\n" + ",\n".join(items) + "\n]")
        else:
            lines.append(head + json.dumps(value, ensure_ascii=False))
    with open(os.path.join(OUT, name), "w", encoding="utf-8") as fh:
        fh.write("{\n" + ",\n".join(lines) + "\n}\n")
    print("wrote", name)


dump("edit_model_cases.json", {
    "how": "MDATK 0.1.1: rename_group, group_references, problems() without selection_problems, describe(); "
           "the project of 'describe' is the input of the rename case 'everything'",
    "rename": RENAME,
    "references": REFERENCES,
    "problems": problems_cases(),
    "describe": describe_cases(),
})
