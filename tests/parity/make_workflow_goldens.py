#!/usr/bin/env python3
"""Project files with planner steps, written by MDATK (tests/data/project-file/workflow).

    python3 -I -B make_workflow_goldens.py MDATK_CHECKOUT OUT_DIR

golden_all_tools: Project() + add_node of every step type + scratch_node of every step type;
golden_edge_cases: unnormalised settings, unknown tools, empty settings, an edge;
golden_test_every_entry / golden_test_saved_in_project: the projects of MDATK's tests
test_projectfile.py::test_every_entry_is_commented and test_workflows.py::test_workflow_saved_in_project;
golden_nondefault (milestone 7, W0-A): one step of every type with every option changed from its default
(check boxes flipped, the last other choice, numbers moved, texts set, the records of repeated sections from
their templates), a scratch record of every type with other settings, files, titles, edges and the naming
keys (use_stamp off, out_dir, step_folders off); golden_written_stamp: the stamp MDATK fixes into the plan
on its first write (Plan.write, a fixed stamp) and the scratch record the wizard keeps.
The "saved" stamp is fixed to 2026-10-15T14:32:00; everything else is what MDATK writes.
"""

import datetime as _dt
import os
import shutil
import sys

if len(sys.argv) != 3:
    sys.exit("usage: make_workflow_goldens.py MDATK_CHECKOUT OUT_DIR")
sys.dont_write_bytecode = True
sys.path.insert(0, sys.argv[1])
OUT = os.path.abspath(sys.argv[2])

from mdatk import projectfile as PF  # noqa: E402
from mdatk.config import PlanEdge, PlanNode, Project  # noqa: E402
from mdatk.workflows import TOOLS  # noqa: E402
from mdatk.workflows import plan as P  # noqa: E402

PF._now = lambda: _dt.datetime(2026, 10, 15, 14, 32, 0)  # noqa: E731


def write(name, project):
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, name), "w", encoding="utf-8", newline="") as fh:
        fh.write(PF.project_text(project))


p = Project()
for k, name in enumerate(TOOLS):
    P.add_node(p.workflow, name, x=40.0 + 280 * k, y=40.0)
for name in TOOLS:
    P.scratch_node(p.workflow, name)
write("golden_all_tools.mdatk.toml", p)

e = Project()
e.workflow.project_name = "edge cases"
n1 = P.add_node(e.workflow, "amber-md", settings={"engine": "pmemd.MPI", "nproc": "4"})
n1.settings["stale_key"] = 1
n1.settings["gpu"] = None
n1.settings["stages"] = n1.settings["stages"][:2]
n1.settings["stages"][1]["stray"] = "x"
n1.files["topology"] = "a.parm7; b.parm7"
n2 = P.add_node(e.workflow, "gromacs")
n2.settings["stages"] = []
e.workflow.nodes.append(PlanNode(id="zz", tool="does-not-exist", step=9, settings={"a": 1, "nested": {"b": 2}}))
e.workflow.nodes.append(PlanNode(id="raw", tool="vmd", step=3, settings={"rna_style": "Tube", "dt": 2}))
e.workflow.nodes.append(PlanNode(id="empty", tool="note", step=0, settings={}))
n3 = P.add_node(e.workflow, "input-file", settings={"path": "maps/x.mrc", "kind": "map"})
n3.title = "Map"
P.connect(e.workflow, n2.id, "topology", P.add_node(e.workflow, "mdatk").id, "topology")
write("golden_edge_cases.mdatk.toml", e)

t = Project()
t.notes = "Aim: compare\nthree \"quoted\" runs"
node = P.add_node(t.workflow, "amber-build", settings={"mg_count": 4})
node.files["structure"] = "model.pdb"
t.layout.panels[1].options = {"smooth": 5, "progressive": True, "custom": 1}
write("golden_test_every_entry.mdatk.toml", t)

s = Project()
s.workflow.project_name = "apt1"
b = P.add_node(s.workflow, "amber-build", settings={"mg_count": 4})
b.files["structure"] = "inputs/model.pdb"
P.add_node(s.workflow, "input-file", settings={"path": "maps/map.mrc", "kind": "map"})
P.connect(s.workflow, b.id, "topology", P.add_node(s.workflow, "amber-md").id, "topology")
P.scratch_node(s.workflow, "vmd").settings["rna_style"] = "Tube"
write("golden_test_saved_in_project.mdatk.toml", s)


def other_value(opt):
    """A value of `opt` other than its default (MDATK's own types, within the option's limits)."""
    if opt.kind == "bool":
        return not opt.default
    if opt.kind in ("radio", "choice"):
        return [v for v, _l in opt.choices if v != opt.default][-1]
    if opt.kind in ("int", "float"):
        if opt.kind == "int":
            value = 7 if opt.default is None else opt.default + 1
        else:
            value = 0.25 if opt.default is None else (opt.default * 1.5 if opt.default else 0.5)
        if opt.maximum is not None and value > opt.maximum:
            value = opt.maximum if opt.default != opt.maximum else opt.maximum - 1
        if opt.minimum is not None and value < opt.minimum:
            value = opt.minimum if opt.default != opt.minimum else opt.minimum + 1
        return int(value) if opt.kind == "int" else float(value)
    if opt.kind == "multiline":
        return f"{opt.key} line 1\n  {opt.key} line 2 \"quoted\""
    if opt.kind == "file":
        return f"inputs/{opt.key}.dat"
    return f"{opt.key} \u00e9 x"


def nondefault(tool):
    out = {}
    for section in tool.sections:
        if section.repeat:
            records = []
            for _tid, _label, values in section.templates:
                rec = section.record_defaults()
                rec.update(values)
                records.append(rec)
            if records:
                rec = dict(records[-1])
                rec.update({o.key: other_value(o) for o in section.options if o.kind != "text"})
                records.append(rec)
            out[section.key] = records
        else:
            out.update({o.key: other_value(o) for o in section.options})
    return out


n = Project()
n.workflow.project_name = "every setting"
n.workflow.use_stamp = False
n.workflow.out_dir = "plans/every"
n.workflow.step_folders = False
for k, name in enumerate(TOOLS):
    node = P.add_node(n.workflow, name, x=40.0 + 280 * (k % 4), y=40.0 + 190 * (k // 4), step=k + 1 if
                      TOOLS[name].generates else 0, title=f"{name} with other settings",
                      settings=nondefault(TOOLS[name]))
    for port in TOOLS[name].input_ports(TOOLS[name].normalise(node.settings)):
        node.files[port.name] = f"inputs/{port.name}.{port.kinds[0]}"
for name in TOOLS:
    P.scratch_node(n.workflow, name).settings.update(nondefault(TOOLS[name]))
ids = {node.tool: node.id for node in n.workflow.nodes}
n.workflow.edges.append(PlanEdge(ids["amber-build"], "topology", ids["amber-md"], "topology"))
n.workflow.edges.append(PlanEdge(ids["farfar2"], "silent", ids["farfar2-models"], "silent"))
write("golden_nondefault.mdatk.toml", n)

w = Project()
w.workflow.project_name = "apt1"
P.apply_template(w.workflow, "amber-gromacs")
w.workflow.nodes[0].files["structure"] = "/golden/inputs/model.pdb"
P.Plan(w.workflow, os.path.join(OUT, "_plan"), w.groups, stamp="20261015-1432").write(None)
shutil.rmtree(os.path.join(OUT, "_plan"))
wz = P.scratch_node(w.workflow, "vmd")
wz.settings["render"] = "snapshot"
wz.files["topology"] = "/golden/inputs/sys.parm7"
write("golden_written_stamp.mdatk.toml", w)
