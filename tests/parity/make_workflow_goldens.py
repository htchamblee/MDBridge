#!/usr/bin/env python3
"""Project files with planner steps, written by MDATK (tests/data/project-file/workflow).

    python3 -I -B make_workflow_goldens.py MDATK_CHECKOUT OUT_DIR

golden_all_tools: Project() + add_node of every step type + scratch_node of every step type;
golden_edge_cases: unnormalised settings, unknown tools, empty settings, an edge;
golden_test_every_entry / golden_test_saved_in_project: the projects of MDATK's tests
test_projectfile.py::test_every_entry_is_commented and test_workflows.py::test_workflow_saved_in_project.
The "saved" stamp is fixed to 2026-10-15T14:32:00; everything else is what MDATK writes.
"""

import datetime as _dt
import os
import sys

if len(sys.argv) != 3:
    sys.exit("usage: make_workflow_goldens.py MDATK_CHECKOUT OUT_DIR")
sys.dont_write_bytecode = True
sys.path.insert(0, sys.argv[1])
OUT = os.path.abspath(sys.argv[2])

from mdatk import projectfile as PF  # noqa: E402
from mdatk.config import PlanNode, Project  # noqa: E402
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
