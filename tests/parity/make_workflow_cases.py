#!/usr/bin/env python3
"""What MDATK v0.0.4.1's workflow layer (the Project Planner's step types, plan, master script and
summary) writes, with MDBridge's improvements (improvements_workflow.py), for tests/workflow
(tests/data/README.md "workflow").

    python3 -I -B tests/parity/make_workflow_cases.py MDATK_CHECKOUT OUT_DIR

OUT_DIR is tests/data/workflow of a tests/data tree (the cases read topology/demo_dry.parm7 next to
it). It writes:

  fixtures/       the input files the plan-time checks read (PDB models built with MDATK's builder,
                  an annotated model, secondary-structure files, a residue library, a frcmod)
  cases/*.json    the cases by group: templates, plan, graph, one file per step type
  texts/<id>.txt  every generated text once, under the first 16 hex digits of its SHA-256 (the
                  cases name their texts by id: a helper written by 40 cases is stored once)

The case design is MDATK 0.1.2's tests/golden.py (the owner's code, ported to 0.0.4.1) with full
texts: every template with dummy inputs; every step type with its defaults, then each check box
flipped and each radio / combo choice taken, one at a time, and the same for the options of the
repeated sections (on one record of the default list, the first whose type the option applies to);
each record template of the repeated sections as a list of its own; segments, production and
continuation stages; typed secondary structures with G-U and U-G pairs (wobble off / on); real input
files where the checks read them; hand-picked paths with spaces (WG-01); plan-level cases (names,
collisions, loops, the master and the summary) and graph edits. All inputs of a step-type case are
connected to dummy absolute paths /golden/inputs/<port>.<kind> (they do not exist); a case with
fixtures copies them into its plan folder (<ROOT>/inputs/...).

Every case runs twice, each time in a child process of this script (--worker): with MDBridge's
improvements (the expected data) and as MDATK v0.0.4.1 (recorded where it differs: "mdatk" in a
case, the issues and the SHA-256 of each text). Texts are masked (improvements_workflow.mask: the
version after the program name as <VERSION>, dates as <DATE>; the plan folder as <ROOT>). The stamp is
20261015-1432. Standard library only in this process; the children import MDATK (numpy for the
fixtures' helices, kernels pinned as in the chem generators). The output is reproducible (also with
other PYTHONHASHSEED values).
"""

from __future__ import annotations

import copy
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)
sys.dont_write_bytecode = True

import chem_common  # noqa: E402

chem_common.pin_kernels()

STAMP = "20261015-1432"
ROOT = "<ROOT>"
DUMMY = "/golden/inputs"
SPACED = "/golden/my inputs"
TYPED_SS = {"sequence": "GGCUCUUCGGGGUC", "dotbracket": "(((((....)))))"}  # G2·U13 and U4·G11
DEMO_SEQUENCE = "GCAUGGCAUCGA UCGAUGCCAUGC"  # topology/demo_dry.parm7: a 12 + 12 duplex (and 3 Mg2+)
DEMO_DOTBRACKET = "(((((((((((( ))))))))))))"
FIXTURE_SEQUENCE = "GGCAUC"


# ================================================================================================
# the cases (built in the worker: the variants follow MDBridge's defaults)

def dummy(port):
    return f"{DUMMY}/{port.name}.{port.kinds[0]}"


def node(nid, tool, step=3, settings=None, files=None, title="", x=40.0, y=40.0):
    return {"id": nid, "tool": tool, "step": step, "title": title, "x": x, "y": y,
            "settings": copy.deepcopy(settings or {}), "files": dict(files or {})}


def plan(nodes, edges=(), name="golden", stamp=STAMP, use_stamp=True, step_folders=True):
    return {"project_name": name, "stamp": stamp, "use_stamp": use_stamp, "step_folders": step_folders,
            "nodes": nodes, "edges": [list(e) for e in edges]}


# the items a step type's default case shows besides D11 (which every generated text shows)
DEFAULT_ITEMS = {"amber-build": ["D3", "D9"], "restraint-dicts": ["D3"], "structure-builder": ["WG-18", "D6"],
                 "farfar2": ["WG-07"], "simrna": ["WG-07", "WG-08", "WG-09"], "farfar2-models": ["WG-09"],
                 "rnamake": ["WG-09"], "gromacs": ["WG-21"]}


def build_cases():
    from mdatk.config import Project
    from mdatk.workflows import TOOLS
    from mdatk.workflows.plan import TEMPLATES, apply_template
    from mdatk.config import WorkflowParams

    groups = [{"name": g.name, "residues": g.residues, "selection": g.selection, "color": g.color}
              for g in Project().groups]
    cases = []

    def add(group, name, workflow, items=(), fixtures=None, groups_=None, mdatk_test="", note=""):
        fixtures = {k: v if "/" in v else f"workflow/fixtures/{v}" for k, v in (fixtures or {}).items()}
        c = {"group": group, "name": name, "items": list(items), "workflow": workflow,
             "groups": groups if groups_ is None else groups_, "fixtures": fixtures}
        if mdatk_test:
            c["mdatk_test"] = mdatk_test
        if note:
            c["note"] = note
        cases.append(c)

    def tool_case(name, label, settings, files=None, optional=True, **kw):
        tool = TOOLS[name]
        s = tool.normalise(copy.deepcopy(settings), [])
        if files is None:
            files = {p.name: dummy(p) for p in tool.input_ports(s) if optional or not p.optional}
        add(name, f"tool:{name}:{label}", plan([node("x", name, 3 if tool.generates else 0, settings, files)]),
            **kw)

    # ---- templates (dummy inputs for every unconnected input, as tests/golden.py) -------------------
    for tname in sorted(TEMPLATES):
        wf = WorkflowParams(project_name="golden", stamp=STAMP)
        apply_template(wf, tname)
        nodes = []
        for n in wf.nodes:
            tool = TOOLS[n.tool]
            files = {}
            for port in tool.input_ports(tool.normalise(n.settings)):
                if not port.optional and not any(e.target == n.id and e.input == port.name for e in wf.edges):
                    files[port.name] = dummy(port)
            nodes.append(node(n.id, n.tool, n.step, _changed(tool, n.settings), files, n.title, n.x, n.y))
        add("templates", f"template:{tname}", plan(nodes, [(e.source, e.output, e.target, e.input) for e in wf.edges]),
            items=["WG-02"] if tname in ("amber", "builder-amber", "farfar2-amber", "amber-gromacs") else [])
    # farfar2-amber with a typed sequence, SimRNA, RNAMake and GROMACS (MDATK's
    # test_all_generated_scripts_pass_bash_syntax_check)
    from mdatk.workflows.plan import add_node, connect
    wf = WorkflowParams(project_name="apt1", stamp=STAMP)
    apply_template(wf, "farfar2-amber")
    for n in wf.nodes:
        if n.tool == "farfar2":
            n.settings.update(sequence="GGGCAC GUGCCC", dotbracket="(((((( ))))))")
    simrna = add_node(wf, "simrna", settings={"sequence": "GGGCAC GUGCCC", "dotbracket": "(((((( ))))))"})
    rm = add_node(wf, "rnamake", settings={"start_bp": "A1-A12", "end_bp": "A6-A7"})
    connect(wf, simrna.id, "best", rm.id, "structure")
    gmx = add_node(wf, "gromacs")
    build = [n for n in wf.nodes if n.tool == "amber-build"][0]
    connect(wf, build.id, "topology", gmx.id, "topology")
    connect(wf, build.id, "coordinates", gmx.id, "coordinates")
    nodes = [node(n.id, n.tool, n.step, _changed(TOOLS[n.tool], n.settings), {}, n.title, n.x, n.y) for n in wf.nodes]
    add("templates", "template:farfar2-amber:bash-syntax-check",
        plan(nodes, [(e.source, e.output, e.target, e.input) for e in wf.edges], name="apt1"),
        mdatk_test="tests/test_workflows.py::test_all_generated_scripts_pass_bash_syntax_check")
    # the amber template written whole (MDATK's test_planner_template_connect_and_write settings)
    wf = WorkflowParams(project_name="apt1", stamp=STAMP)
    apply_template(wf, "amber")
    nodes = []
    for n in wf.nodes:
        tool = TOOLS[n.tool]
        settings = _changed(tool, n.settings)
        files = {}
        if n.tool == "amber-build":
            files = {"structure": f"{ROOT}/inputs/model.pdb"}
        if n.tool == "amber-restraints":
            settings.update(sequence=FIXTURE_SEQUENCE + " GAUGCC", dotbracket="(((((( ))))))")
        nodes.append(node(n.id, n.tool, n.step, settings, files, n.title, n.x, n.y))
    add("templates", "template:amber:written", plan(nodes, [(e.source, e.output, e.target, e.input) for e in wf.edges],
                                                    name="apt1"),
        fixtures={"inputs/model.pdb": "d_duplex.pdb"}, items=["WG-06"],
        mdatk_test="tests/test_workflow_gui.py::test_planner_template_connect_and_write")

    # ---- every step type: defaults, then one change at a time ------------------------------------
    for name in sorted(TOOLS):
        tool = TOOLS[name]
        for label, change in variants(tool):
            items = DEFAULT_ITEMS.get(name, []) if label == "default" else []
            items += {("vmd", "render=movie"): ["WG-09"], ("amber-build", "split_types=False"): ["D9"],
                      ("simrna", "cluster=False"): ["WG-20"], ("farfar2-models", "cluster=False"): ["WG-09", "WG-19"],
                      ("farfar2-models", "default"): ["WG-19"]}.get((name, label), [])
            tool_case(name, label, change, items=items)
        for label, change in record_variants(tool):
            tool_case(name, label, change)
        if name == "vmd":
            add(name, "tool:vmd:default:no-groups",
                plan([node("x", "vmd", 3, {}, {p.name: dummy(p) for p in tool.input_ports(tool.normalise({}))})]),
                groups_=[])
    extra_tool_cases(tool_case, add, TOOLS)
    plan_cases(add, TOOLS)
    return cases


def _changed(tool, settings):
    """The settings that differ from the defaults (what a case gives; the plan normalises)."""
    defaults = tool.normalise({})
    s = tool.normalise(copy.deepcopy(settings), [])
    return {k: v for k, v in s.items() if v != defaults.get(k)}


def variants(tool):
    """(label, settings): the defaults, then one changed value at a time (tests/golden.py _variants)."""
    yield "default", {}
    for section in tool.sections:
        if section.repeat:
            continue
        for opt in section.options:
            if opt.kind == "bool":
                yield f"{opt.key}={not opt.default}", {opt.key: not opt.default}
            elif opt.kind in ("radio", "choice"):
                for value, _label in opt.choices:
                    if value != opt.default:
                        yield f"{opt.key}={value}", {opt.key: value}


def record_variants(tool):
    """The repeated sections: each record template as a list of its own, then each check box / choice of
    the record options changed on one record of the default list (the first the option applies to)."""
    from mdatk.workflows.schema import requirements_met

    for section in tool.sections:
        if not section.repeat:
            continue
        defaults = section.default_records()
        for tid, _label, values in section.templates:
            rec = section.record_defaults()
            rec.update(copy.deepcopy(values))
            yield f"{section.key}=[{tid}]", {section.key: [rec]}
        for opt in section.options:
            if opt.kind == "bool":
                changes = [not opt.default]
            elif opt.kind in ("radio", "choice"):
                changes = [v for v, _l in opt.choices]
            else:
                continue
            for value in changes:
                k = next((k for k, rec in enumerate(defaults)
                          if rec.get(opt.key) != value and requirements_met(opt.requires, rec)), None)
                if k is None:
                    continue
                records = copy.deepcopy(defaults)
                records[k][opt.key] = value
                yield f"{section.key}.{k}.{opt.key}={value}", {section.key: records}


def extra_tool_cases(tool_case, add, TOOLS):
    """Typed secondary structures, segments and production, fixtures, numbers, spaced paths."""
    T = TOOLS
    ss = dict(TYPED_SS)
    for name, more in (("amber-restraints", {}), ("phenix-rsr", {"resolution": 3.5}), ("simrna", {}),
                       ("farfar2", {})):
        if name in ("amber-restraints", "phenix-rsr"):
            tool_case(name, "typed-ss", dict(ss, wobble=False, **more), optional=False)
            tool_case(name, "typed-ss+wobble", dict(ss, wobble=True, **more), optional=False)
        else:
            tool_case(name, "typed-ss", dict(ss, **more), optional=False)
    tool_case("structure-builder", "mode=fit+resolution", {"mode": "fit", "resolution": 4.0})
    tool_case("farfar2", "typed-ss+stems", dict(ss, input_res_mode="stems"), optional=False)
    tool_case("farfar2-models", "n_models=105", {"n_models": 105, "cluster": False}, items=["WG-19"])
    tool_case("farfar2-models", "n_models=105+json", {"n_models": 105, "dssr_json": True}, items=["WG-19"])

    # amber-md: segments, production, continuation, ramps (WG-04, WG-05)
    md = T["amber-md"]
    base = md.normalise({})["stages"]

    def stages(*edits):
        records = copy.deepcopy(base)
        for k, values in edits:
            records[k].update(values)
        return records

    tool_case("amber-md", "segments:heat=3", {"stages": stages((2, {"segments": 3}))}, items=["WG-04"])
    tool_case("amber-md", "segments:heat=3+ramp", {"stages": stages((2, {"segments": 3, "ramp": True,
                                                                         "ramp_steps": 20000}))}, items=["WG-04"])
    tool_case("amber-md", "segments:eq=2+ramp-continued",
              {"stages": stages((3, {"segments": 2, "ramp": True, "tempi": 280.0, "ramp_steps": 10000}))},
              items=["WG-04"])
    tool_case("amber-md", "segments:prod=4", {"stages": stages((5, {"segments": 4}))})
    tool_case("amber-md", "segments:eq=2+stage-eq-001",
              {"stages": stages((3, {"segments": 2}), (4, {"name": "eq-001"}))}, items=["WG-05"])
    tool_case("amber-md", "segments:heat=2+stage-heat-cont",
              {"stages": stages((2, {"segments": 2}), (4, {"name": "heat-cont"}))}, items=["WG-05"])
    tool_case("amber-md", "segments:eq=2+eq-001=2",
              {"stages": stages((3, {"segments": 2}), (4, {"name": "eq-001", "segments": 2}))}, items=["WG-05"])
    tool_case("amber-md", "two-production-stages",
              {"stages": stages((4, {"production": True}))})
    tool_case("amber-md", "no-production+separate", {"stages": stages((5, {"production": False}))})
    tool_case("amber-md", "continue-after-min", {"stages": stages((2, {"restart": "continue"}))})
    tool_case("amber-md", "continue-first+input-velocities",
              {"input_velocities": True, "stages": [dict(base[2], restart="continue")]})
    tool_case("amber-md", "dt=0.004-no-hmr", {"stages": stages((3, {"dt": 0.004}))})
    tool_case("amber-md", "dt=0.004-no-shake", {"stages": stages((3, {"dt": 0.004, "shake": False}))})
    tool_case("amber-md", "numbers", {"engine": "pmemd.MPI", "nproc": 16, "mpi_launcher": "srun -n {n}",
                                      "stages": stages((0, {"maxcyc": 5000, "ncyc": 6000}),
                                                       (3, {"tol": 1e-7, "ig": 12345, "cut": 9.0, "pres0": 1.01325,
                                                            "segments": 2, "extra_cntrl": "nscm=1000\nntxo=2,",
                                                            "extra_namelists": "&ewald\n  skinnb=2.0,\n/"}),
                                                       (5, {"nstlim": 1234567, "ntwx": 1000}))})
    tool_case("amber-md", "stages=[]", {"stages": []})
    tool_case("amber-md", "posres-mask", {"stages": stages((0, {"posres_mask": ":1-24 & !@H="}))})
    tool_case("amber-md", "gpu+env", {"gpu": "1", "conda_env": "AmberTools26", "env_setup": "module load cuda\n"})
    for value in ("nan", "inf", "-inf"):
        tool_case("amber-build", f"buffer={value}", {"buffer": value}, items=["WG-11", "D13"])
        tool_case("amber-md", f"stages.3.dt={value}", {"stages": stages((3, {"dt": value}))}, items=["WG-11", "D13"])
    tool_case("amber-md", "stages.0.maxcyc=inf", {"stages": stages((0, {"maxcyc": "inf"}))}, items=["WG-12"])
    tool_case("amber-md", "stages.0.maxcyc=nan", {"stages": stages((0, {"maxcyc": "nan"}))}, items=["WG-12"])

    # gromacs: numbers
    g = T["gromacs"]
    gst = g.normalise({})["stages"]
    tool_case("gromacs", "numbers", {"stages": [dict(gst[0], nsteps=12345)] + copy.deepcopy(gst[1:])})

    # amber-build: options with values
    tool_case("amber-build", "ions", {"mg_count": 20, "salt_pairs": 15, "buffer": 12.5, "cation": "K+"})
    for model in ("micro", "nano"):
        for water in ("tip3p", "opc", "tip4p2005"):
            tool_case("amber-build", f"mg_model={model}+water={water}+mg_count=4",
                      {"mg_model": model, "water": water, "mg_count": 4}, items=["D9"])
        tool_case("amber-build", f"mg_model={model}+two_pass+lrna_types",
                  {"mg_model": model, "mg_count": 4, "two_pass": True, "lrna_types": "LD LH LE LX:OH",
                   "parmed_extra": "changeLJPair :MG @%LD 2.0 0.1\nchangeLJPair @%XA @%YB 3.0 0.2"}, items=["D9"])
        tool_case("amber-build", f"mg_model={model}+lrna_types-missing-ref",
                  {"mg_model": model, "mg_count": 4, "lrna_types": "LD LQ"}, items=["D9"])
        tool_case("amber-build", f"mg_model={model}+extra_libs",
                  {"mg_model": model, "mg_count": 4, "extra_libs": "LRNA.lib", "extra_frcmod": "frcmod.LRNA",
                   "cp_source": "frcmod"}, items=["D9"])
    tool_case("amber-build", "lrna_types=bad", {"lrna_types": "LD toolong L:X:Y"}, items=["D9"])
    tool_case("amber-build", "lrna=False+lrna_types", {"lrna": False, "nbfix19": True}, items=["D9"])
    tool_case("amber-build", "two_pass+parmed_extra-mg", {"two_pass": True, "mg_count": 2,
                                                          "parmed_extra": "changeLJPair :MG @%OS 2.1 0.2"},
              items=["D9"])
    tool_case("amber-build", "leap-lines", {"extra_frcmod": "frcmod.LRNA\nfrcmod.extra", "extra_libs": "LRNA.lib",
                                            "leap_before_params": 'addAtomTypes { { "LD" "O" "sp3" } }',
                                            "leap_before_load": 'addPdbResMap { { "LG" "G0" } }',
                                            "leap_after_load": "bond mol.1.P mol.2.O3'", "lrna_ranges": "B"})
    tool_case("amber-build", "lrna_ranges=bad", {"lrna_ranges": "A:4-1"})
    tool_case("amber-build", "env", {"conda_env": "AmberTools26", "env_setup": "source /opt/amber26/amber.sh"})

    # plan-time reads: fixtures in the plan folder
    pdb = {"structure": f"{ROOT}/inputs/model.pdb"}
    for fixture, items in (("d_duplex.pdb", []), ("l_duplex.pdb", ["D9"]), ("dl_duplex.pdb", ["D9"]),
                           ("blank_chains.pdb", ["WG-06"])):
        stem = fixture[:-4]
        fx = {"inputs/model.pdb": fixture}
        tool_case("amber-build", f"fixture:{stem}", {}, files=pdb, fixtures=fx, items=items)
        tool_case("amber-build", f"fixture:{stem}:libs", {"extra_libs": f"{ROOT}/inputs/LRNA.lib",
                                                         "extra_frcmod": "frcmod.LRNA"},
                  files=pdb, fixtures=dict(fx, **{"inputs/LRNA.lib": "lrna_a0.lib"}), items=items)
    for label, settings in (
            ("libs+resmap", {"leap_before_load": 'addPdbResMap { { "LG5" "G05" } { "LC3" "C03" } }\n'
                                                 'addPdbResMap {\n  { "LG" "G0" } { "LC" "C0" }\n  { "LA" "A0" } '
                                                 '{ "LU" "U0" }\n}'}),
            ("libs+resmap-partial", {"leap_before_load": 'addPdbResMap { { "LG" "G0" } }'}),
            ("libs+assigned", {"leap_before_params": "LG = copy G0\nLC = copy C0\nLG5 = copy G05\nLC3 = copy C03"}),
            ("libs+loadoff", {"extra_libs": "", "leap_before_load": f'loadOff "{ROOT}/inputs/LRNA.lib"'}),
            ("libs-unread", {"extra_libs": "nucleic12_LRNA.lib"}),
            ("libs+no-terminals", {"lrna_terminals": False}),
            ("libs+lrna_ranges", {"lrna_ranges": "A"}),
            ("libs+mg_model=nano", {"mg_model": "nano", "mg_count": 2})):
        s = {"extra_libs": f"{ROOT}/inputs/LRNA.lib", "extra_frcmod": "frcmod.LRNA"}
        s.update(settings)
        tool_case("amber-build", f"fixture:l_duplex:{label}", s, files=pdb,
                  fixtures={"inputs/model.pdb": "l_duplex.pdb", "inputs/LRNA.lib": "lrna_a0.lib"}, items=["D9"])
    tool_case("amber-build", "fixture:l_duplex:libs-dudek", {"extra_libs": f"{ROOT}/inputs/LRNA.lib",
                                                             "extra_frcmod": "frcmod.LRNA"}, files=pdb,
              fixtures={"inputs/model.pdb": "l_duplex.pdb", "inputs/LRNA.lib": "lrna_dudek.lib"}, items=["D9"])
    tool_case("amber-build", "fixture:cp_frcmod-new-types", {"cp_source": "frcmod",
                                                             "cp_frcmod": f"{ROOT}/inputs/frcmod.cp"},
              files=pdb, fixtures={"inputs/model.pdb": "d_duplex.pdb", "inputs/frcmod.cp": "cp_types.frcmod"})
    # Amber build -> MD: the solute count of the input PDB in the default restraint mask (WG-06)
    for fixture in ("d_duplex.pdb", "blank_chains.pdb"):
        add("amber-md", f"plan:amber-build+md:{fixture[:-4]}",
            plan([node("b", "amber-build", 1, {}, pdb), node("m", "amber-md", 2, {}, {"restraints": f"{DUMMY}/h.dist"},
                                                             x=320.0)],
                 [("b", "topology", "m", "topology"), ("b", "coordinates", "m", "coordinates")]),
            fixtures={"inputs/model.pdb": fixture}, items=["WG-06"])

    # amber-restraints: the direct DISANG file of an existing topology, secondary-structure files
    top = {"topology": f"{ROOT}/inputs/demo_dry.parm7"}
    demo = {"inputs/demo_dry.parm7": "topology/demo_dry.parm7"}
    tool_case("amber-restraints", "fixture:demo-topology",
              {"sequence": DEMO_SEQUENCE, "dotbracket": DEMO_DOTBRACKET}, files=top, fixtures=demo)
    tool_case("amber-restraints", "fixture:demo-topology+offset",
              {"sequence": DEMO_SEQUENCE, "dotbracket": DEMO_DOTBRACKET, "res_offset": 1, "wobble": True},
              files=top, fixtures=demo)
    tool_case("amber-restraints", "fixture:demo-topology-wrong-sequence",
              {"sequence": "GGGGGGGGGGGG CCCCCCCCCCCC", "dotbracket": DEMO_DOTBRACKET}, files=top, fixtures=demo)
    tool_case("amber-restraints", "fixture:demo-topology-window",
              {"sequence": DEMO_SEQUENCE, "dotbracket": DEMO_DOTBRACKET, "r1": 1.5, "r2": 2.25, "r3": 3.125,
               "r4": 4.0, "rk2": 20.0, "rk3": 7.5}, files=top, fixtures=demo)
    tool_case("amber-restraints", "r-order", dict(TYPED_SS, r1=3.0, r2=2.0), optional=False)
    for name, port in (("amber-restraints", "secstruct"), ("simrna", "secstruct"), ("phenix-rsr", "secstruct")):
        tool = T[name]
        base_files = {p.name: dummy(p) for p in tool.input_ports(tool.normalise({})) if not p.optional}
        more = {"resolution": 3.5} if name == "phenix-rsr" else {}
        for fixture, items in (("hairpin.secstruct", []), ("hairpin_no_sequence.secstruct", ["WG-10"])):
            files = dict(base_files, **{port: f"{ROOT}/inputs/{fixture}"})
            tool_case(name, f"fixture:{fixture[:-10]}", dict(more), files=files,
                      fixtures={f"inputs/{fixture}": fixture}, items=items)
        files = dict(base_files, **{port: f"{ROOT}/inputs/hairpin_no_sequence.secstruct"})
        tool_case(name, "fixture:hairpin_no_sequence+typed", dict(more, sequence=TYPED_SS["sequence"]), files=files,
                  fixtures={"inputs/hairpin_no_sequence.secstruct": "hairpin_no_sequence.secstruct"},
                  items=["WG-10"])
    # restraint dictionaries: the model's codes
    for fixture in ("d_duplex.pdb", "l_duplex.pdb", "dl_duplex.pdb"):
        tool_case("restraint-dicts", f"fixture:{fixture[:-4]}", {}, files={"structure": f"{ROOT}/inputs/model.pdb"},
                  fixtures={"inputs/model.pdb": fixture})
    tool_case("restraint-dicts", "custom-codes", {"codes_mode": "list", "codes": "0G5 0A XYZ", "custom":
                                                  "XYZ CCO"}, optional=False)
    # phenix-rsr / coot-rsr: the model's annotation and names
    for name in ("phenix-rsr", "coot-rsr"):
        tool = T[name]
        base_files = {p.name: dummy(p) for p in tool.input_ports(tool.normalise({})) if not p.optional}
        for fixture in ("annotated.pdb", "l_duplex.pdb", "d_duplex.pdb"):
            more = {"resolution": 3.5} if name == "phenix-rsr" else {}
            tool_case(name, f"fixture:{fixture[:-4]}", more,
                      files=dict(base_files, model=f"{ROOT}/inputs/model.pdb"), fixtures={"inputs/model.pdb": fixture})
    # structure builder: mutations against the model
    tool_case("structure-builder", "fixture:mutate", {"mode": "mutate", "mutations": "A:2C, B:12G"},
              files={"pdb": f"{ROOT}/inputs/model.pdb"}, fixtures={"inputs/model.pdb": "d_duplex.pdb"})
    tool_case("structure-builder", "fixture:mutate-missing", {"mode": "mutate", "mutations": "A:2C, C:7G"},
              files={"pdb": f"{ROOT}/inputs/model.pdb"}, fixtures={"inputs/model.pdb": "d_duplex.pdb"})
    tool_case("structure-builder", "python=python3.12", {"python": "python3.12", "sequence": "GGCAUC"},
              items=["WG-18"])
    tool_case("structure-builder", "helix-values", {"sequence": "GGC AUC", "chains": "X,Y", "start1": -5,
                                                    "start2": 101, "scheme": "ccd"}, items=["WG-18"])
    # input-file: an existing file
    tool_case("input-file", "fixture:exists", {"path": f"{ROOT}/inputs/model.pdb", "kind": "pdb"}, files={},
              fixtures={"inputs/model.pdb": "d_duplex.pdb"})
    tool_case("input-file", "missing", {"path": "/golden/inputs/none.map", "kind": "map"}, files={})
    tool_case("note", "text", {"text": "check the Mg count\n\nthen run"}, files={})

    # WG-01: hand-picked paths with spaces
    for name in sorted(T):
        tool = T[name]
        if not tool.generates:
            continue
        s = tool.normalise({})
        files = {p.name: f"{SPACED}/{p.name}.{p.kinds[0]}" for p in tool.input_ports(s)}
        if not files:
            continue
        settings = {}
        if name == "amber-md":
            settings = {"restraint_ref": "initial"}
        if name == "farfar2-models":
            files["silent"] = f"{SPACED}/a.out;{SPACED}/b c.out"
        tool_case(name, "spaced-paths", settings, files=files, items=["WG-01"])
    tool_case("amber-md", "spaced-restraints-only", {}, files={"topology": f"{DUMMY}/topology.parm7",
                                                               "coordinates": f"{DUMMY}/coordinates.rst7",
                                                               "restraints": f"{SPACED}/helix.dist"}, items=["WG-01"])
    tool_case("amber-md", "quote-in-path", {}, files={"topology": f"{DUMMY}/topology.parm7",
                                                      "coordinates": "/golden/it's/sys.rst7",
                                                      "restraints": f"{DUMMY}/h.dist"}, items=["WG-01"])


def plan_cases(add, TOOLS):
    """Plan-level behaviour: naming, collisions, loops, edges, the master and the summary."""
    d = lambda port, kind: f"{DUMMY}/{port}.{kind}"  # noqa: E731
    build_files = {"structure": d("structure", "pdb")}
    md_files = {"restraints": d("h", "dist")}
    ss = dict(TYPED_SS)

    add("plan", "plan:empty", plan([]))
    add("plan", "plan:no-stamp", plan([node("b", "amber-build", 1, {}, build_files)], use_stamp=False))
    add("plan", "plan:no-step-folders", plan([node("b", "amber-build", 1, {}, build_files)], step_folders=False))
    add("plan", "plan:project-name", plan([node("b", "amber-build", 1, {}, build_files)], name="  my plan!/v2. "))
    add("plan", "plan:project-name-empty", plan([node("n", "note", 0, {"text": "x"})], name="..."))
    add("plan", "plan:unknown-tool", plan([node("a", "does-not-exist", 1), node("b", "amber-build", 2, {}, build_files)]))
    add("plan", "plan:duplicate-id", plan([node("a", "amber-build", 1, {}, build_files),
                                           node("a", "amber-md", 2, {}, md_files), node("", "note", 0)]))
    add("plan", "plan:loop", plan([node("a", "amber-md", 1, {"stages": []}), node("b", "amber-md", 2)],
                                  [("a", "topology", "b", "topology"), ("b", "topology", "a", "topology")]))
    add("plan", "plan:missing-step-edge", plan([node("a", "amber-build", 1, {}, build_files)],
                                               [("a", "topology", "zz", "topology")]))
    add("plan", "plan:port-no-longer-produced",
        plan([node("a", "amber-build", 1, {"save_pdb": False}, build_files), node("v", "vmd", 2)],
             [("a", "pdb", "v", "topology")]))
    add("plan", "plan:input-no-longer-exists",
        plan([node("a", "amber-build", 1, {}, build_files), node("v", "vmd", 2, {}, {"topology": d("t", "parm7")})],
             [("a", "topology", "v", "nothing")]))
    add("plan", "plan:kind-mismatch",
        plan([node("a", "amber-build", 1, {}, build_files), node("m", "amber-md", 2, {}, md_files)],
             [("a", "pdb", "m", "topology"), ("a", "coordinates", "m", "coordinates")]))
    add("plan", "plan:two-into-single",
        plan([node("a", "amber-build", 1, {}, build_files), node("b", "amber-build", 2, {}, build_files),
              node("m", "amber-md", 3, {}, md_files)],
             [("a", "topology", "m", "topology"), ("b", "topology", "m", "topology"),
              ("a", "coordinates", "m", "coordinates")]))
    add("plan", "plan:connected-and-chosen",
        plan([node("a", "amber-build", 1, {}, build_files),
              node("m", "amber-md", 2, {}, dict(md_files, topology=d("other", "parm7")))],
             [("a", "topology", "m", "topology"), ("a", "coordinates", "m", "coordinates")]))
    add("plan", "plan:not-produced-with-settings",
        plan([node("m", "amber-md", 1, {"stages": [dict(TOOLS["amber-md"].normalise({})["stages"][0])]},
                   {"topology": d("t", "parm7"), "coordinates": d("c", "rst7")}),
              node("v", "vmd", 2)],
             [("m", "topology", "v", "topology"), ("m", "trajectory", "v", "trajectory")]))
    add("plan", "plan:multiple-into-multiple",
        plan([node("f1", "farfar2", 1, dict(ss)), node("f2", "farfar2", 2, dict(ss)),
              node("m", "farfar2-models", 3)],
             [("f1", "silent", "m", "silent"), ("f2", "silent", "m", "silent")]))
    add("plan", "plan:step-zero", plan([node("a", "amber-build", 0, {}, build_files)]))
    add("plan", "plan:same-step-and-tag", plan([node("a", "amber-build", 2, {}, build_files),
                                                node("b", "amber-build", 2, {}, build_files)]))
    add("plan", "plan:shared-number", plan([node("a", "amber-build", 2, {}, build_files),
                                            node("r", "amber-restraints", 2, dict(ss), {"topology": d("t", "parm7")})]))
    # WG-03: SimRNA and FARFAR2 models at step 02 without step folders promise the same model files
    shared = [node("s", "simrna", 2, dict(ss)), node("m", "farfar2-models", 2, {}, {"silent": d("a", "out")})]
    add("plan", "plan:shared-outputs", plan(copy.deepcopy(shared), step_folders=False), items=["WG-03"])
    add("plan", "plan:shared-outputs+step-folders", plan(copy.deepcopy(shared)), items=["WG-03"])
    add("plan", "plan:order-by-position",
        plan([node("c", "note", 0, {"text": "c"}, x=500.0, y=40.0), node("b", "note", 0, {"text": "b"}, x=40.0,
                                                                          y=300.0),
              node("a", "note", 0, {"text": "a"}, x=300.0, y=40.0),
              node("r2", "amber-restraints", 1, dict(ss), {"topology": d("t", "parm7")}, x=40.0, y=500.0),
              node("r1", "amber-restraints", 1, dict(ss), {"topology": d("u", "parm7")}, x=900.0, y=10.0)]))
    # WG-02: a step with errors, the steps reading from it, an independent one
    add("plan", "plan:errors-upstream",
        plan([node("b", "amber-build", 1, {}, {}), node("r", "amber-restraints", 2, dict(ss)),
              node("m", "amber-md", 3), node("v", "vmd", 4), node("i", "amber-restraints", 5, dict(ss),
                                                                    {"topology": d("t", "parm7")})],
             [("b", "topology", "r", "topology"), ("b", "topology", "m", "topology"),
              ("b", "coordinates", "m", "coordinates"), ("r", "restraints", "m", "restraints"),
              ("m", "topology", "v", "topology"), ("m", "trajectory", "v", "trajectory")]), items=["WG-02"])
    add("plan", "plan:errors-two-sources",
        plan([node("a", "amber-build", 1, {}, {}), node("r", "amber-restraints", 2, {}, {"topology": d("t", "parm7")}),
              node("m", "amber-md", 3)],
             [("a", "topology", "m", "topology"), ("a", "coordinates", "m", "coordinates"),
              ("r", "restraints", "m", "restraints")]), items=["WG-02"])
    add("plan", "plan:notes-and-analysis",
        plan([node("n", "note", 0, {"text": "Aim: compare\n\nthree \"quoted\" runs"}),
              node("f", "input-file", 0, {"path": "/golden/inputs/x.mrc", "kind": "map"}),
              node("a", "amber-build", 1, {}, build_files), node("m", "amber-md", 2, {}, md_files),
              node("k", "mdatk", 0)],
             [("a", "topology", "m", "topology"), ("a", "coordinates", "m", "coordinates"),
              ("m", "topology", "k", "topology"), ("m", "trajectory", "k", "trajectory")]))
    add("plan", "plan:titles", plan([node("a", "amber-build", 1, {}, build_files, title="System build (L-RNA)"),
                                     node("n", "note", 0, {"text": "x"}, title="Remember")]))
    add("plan", "plan:mdbridge-steps",
        plan([node("h", "structure-builder", 1, {"sequence": "GGCAUC"}), node("b", "amber-build", 2),
              node("d", "restraint-dicts", 3, {}, {"structure": d("model", "pdb")})],
             [("h", "pdb", "b", "structure")]), items=["D3", "D6", "WG-18"])


# ================================================================================================
# the worker: run the cases with MDATK (improved or raw)

def jsonable(value):
    import dataclasses

    if dataclasses.is_dataclass(value) and not isinstance(value, type):
        return {f.name: jsonable(getattr(value, f.name)) for f in dataclasses.fields(value)}
    if isinstance(value, dict):
        return {str(k): jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, set, frozenset)):
        items = [jsonable(v) for v in value]
        return sorted(items, key=json.dumps) if isinstance(value, (set, frozenset)) else items
    if isinstance(value, float) and (value != value or value in (float("inf"), float("-inf"))):
        return repr(value)
    return value


class Texts:
    def __init__(self):
        self.store = {}

    def ref(self, text):
        key = hashlib.sha256(text.encode("utf-8", "surrogateescape")).hexdigest()[:16]
        old = self.store.setdefault(key, text)
        if old != text:
            raise RuntimeError(f"text id collision {key}")
        return key


def worker(kind, mdatk, cases_file, out_file, data_dir):
    import logging

    logging.disable(logging.CRITICAL)  # normalise() logs the unknown keys and bad values it drops
    sys.path.insert(0, os.path.abspath(mdatk))
    import improvements_workflow as IW

    if kind == "improved":
        IW.apply_mdatk()
    from mdatk.config import PlanEdge, PlanNode, ResidueGroup, WorkflowParams
    from mdatk.workflows import TOOLS
    from mdatk.workflows import plan as P
    from mdatk.workflows.schema import Issue

    if cases_file is None:
        cases = build_cases() + graph_cases()
    else:
        with open(cases_file, encoding="utf-8") as fh:
            cases = json.load(fh)
    texts = Texts()
    results = []
    for c in cases:
        if c["group"] == "graph":
            results.append(run_graph(c, TOOLS, P))
            continue
        tmp = tempfile.mkdtemp(prefix="mdbridge-workflow-")
        try:
            root = os.path.join(tmp, "plan")
            os.makedirs(root)
            for dest, source in c["fixtures"].items():
                target = os.path.join(root, dest)
                os.makedirs(os.path.dirname(target), exist_ok=True)
                shutil.copyfile(os.path.join(data_dir, source), target)
            try:
                results.append(run_case(c, root, texts, IW, PlanNode, PlanEdge, WorkflowParams, ResidueGroup, P,
                                        Issue))
            except Exception as exc:  # noqa: BLE001 - MDATK crashes on some (WG-11): recorded
                if kind == "improved":
                    raise
                results.append({"name": c["name"], "crash": f"{type(exc).__name__}: {exc}"})
        finally:
            shutil.rmtree(tmp, ignore_errors=True)
    with open(out_file, "w", encoding="utf-8") as fh:
        json.dump({"cases": cases, "results": results, "texts": texts.store}, fh, ensure_ascii=False)


def run_case(c, root, texts, IW, PlanNode, PlanEdge, WorkflowParams, ResidueGroup, P, Issue):
    def unroot(text):
        return IW.mask(text.replace(root, ROOT))

    def deroot(value):
        if isinstance(value, str):
            return value.replace(ROOT, root)
        if isinstance(value, list):
            return [deroot(v) for v in value]
        if isinstance(value, dict):
            return {k: deroot(v) for k, v in value.items()}
        return value

    def back(value):
        if isinstance(value, str):
            return unroot(value)
        if isinstance(value, list):
            return [back(v) for v in value]
        if isinstance(value, dict):
            return {k: back(v) for k, v in value.items()}
        return value

    w = c["workflow"]
    wf = WorkflowParams(project_name=w["project_name"], stamp=w["stamp"], use_stamp=w["use_stamp"],
                        step_folders=w["step_folders"],
                        nodes=[PlanNode(id=n["id"], tool=n["tool"], step=n["step"], title=n["title"], x=n["x"],
                                        y=n["y"], settings=deroot(copy.deepcopy(n["settings"])),
                                        files=deroot(dict(n["files"]))) for n in w["nodes"]],
                        edges=[PlanEdge(*e) for e in w["edges"]])
    groups = [ResidueGroup(g["name"], g["residues"], g["selection"], g["color"]) for g in c["groups"]]
    plan = P.Plan(wf, root, groups)

    def issue(i):
        return [i.level, unroot(i.message), list(i.keys)]

    def port_value(pv):
        return back({"paths": list(pv.paths), "kind": pv.kind, "data": jsonable(pv.data), "source": pv.source,
                     "source_label": pv.source_label})

    out = {"name": c["name"], "issues": [[nid, *issue(i)] for nid, i in plan.all_issues()], "order": list(plan.order),
           "steps": {}}
    for nid in plan.order + [n for n in plan.steps if n not in plan.order]:
        step = plan.steps[nid]
        rec = {"label": step.label, "ok": step.ok}
        if step.tool is not None:
            rec["describe"] = unroot(step.tool.describe(step.settings))
            rec["resources"] = list(step.tool.resources(step.settings))
            if step.ctx is not None:
                rec["folder"] = step.ctx.folder
                rec["inputs"] = {k: port_value(v) for k, v in step.ctx.inputs.items()}
            rec["outputs"] = {k: port_value(v) for k, v in step.outputs.items()}
        out["steps"][nid] = rec
    report = plan.write(None)
    for nid in plan.order:
        step = plan.steps[nid]
        if step.tool is None or not step.tool.generates or step.ctx is None:
            continue
        files = step.files
        rec = out["steps"][nid]
        if not files and not step.ok:  # the generator's files although the step has errors (tests/golden.py)
            try:
                files = P.generate_step(step)
                rec["forced"] = True
            except Exception as exc:  # noqa: BLE001 - recorded
                rec["forced"] = True
                rec["exception"] = type(exc).__name__
                files = []
        rec["files"] = []
        for f in files:
            entry = {"name": f.name, "role": f.role, "description": unroot(f.description),
                     "text": texts.ref(unroot(f.content if f.content.endswith("\n") else f.content + "\n"))}
            if not f.content.endswith("\n"):
                entry["eol_added"] = True
            rec["files"].append(entry)
    written = {}
    for path in report.written:
        rel = os.path.relpath(path, root).replace(os.sep, "/")
        with open(path, encoding="utf-8", errors="surrogateescape", newline="") as fh:
            data = fh.read()
        written[rel] = {"text": texts.ref(unroot(data)), "mode": oct(os.stat(path).st_mode & 0o777)[2:]}
    out["write"] = {"written": written, "skipped": {k: unroot(v) for k, v in report.skipped.items()},
                    "master": os.path.relpath(report.master, root) if report.master else "",
                    "summary": os.path.relpath(report.summary, root) if report.summary else "",
                    "stamp": wf.stamp}
    return out


# ---- graph edits ---------------------------------------------------------------------------------

def graph_cases():
    """Graph edits of the planner (plan.add_node, connect, remove_node, apply_template, new_node_id,
    next_step, scratch_node, single_step_plan, analysis_inputs) as operation lists; each operation's
    result or error."""
    ops = []

    def add(name, steps, items=()):
        ops.append({"group": "graph", "name": name, "items": list(items), "ops": steps, "workflow": None,
                    "groups": [], "fixtures": {}})

    add("graph:add-and-connect", [
        ["add_node", "amber-build", 40, 40], ["add_node", "amber-restraints", 320, 230],
        ["add_node", "amber-md", 600, 40], ["add_node", "note", 40, 400], ["add_node", "vmd", 880, 40],
        ["connect", "n1", "topology", "n2", "topology"], ["connect", "n1", "topology", "n3", "topology"],
        ["connect", "n1", "coordinates", "n3", "coordinates"], ["connect", "n2", "restraints", "n3", "restraints"],
        ["connect", "n1", "pdb", "n3", "topology"], ["connect", "n3", "topology", "n1", "structure"],
        ["connect", "n3", "trajectory", "n3", "coordinates"], ["connect", "n1", "nothing", "n3", "topology"],
        ["connect", "n1", "topology", "zz", "topology"], ["connect", "n1", "pdb", "n5", "topology"],
        ["connect", "n3", "topology", "n5", "topology"], ["connect", "n3", "trajectory", "n5", "trajectory"],
        ["connect", "n3", "trajectory", "n5", "trajectory"], ["next_step"], ["new_node_id"],
        ["remove_node", "n2"], ["new_node_id"], ["add_node", "farfar2", 40, 600], ["add_node", "unknown", 0, 0],
    ])
    add("graph:multiple-input", [
        ["add_node", "farfar2", 40, 40], ["add_node", "farfar2", 40, 230], ["add_node", "farfar2-models", 320, 40],
        ["connect", "n1", "silent", "n3", "silent"], ["connect", "n2", "silent", "n3", "silent"],
        ["connect", "n1", "silent", "n3", "silent"]])
    add("graph:templates", [["apply_template", "amber", 40, 40], ["apply_template", "amber-gromacs", 40, 500],
                            ["apply_template", "cryoem-refine", 40, 900], ["apply_template", "nothing", 0, 0],
                            ["next_step"]])
    add("graph:scratch", [["scratch_node", "amber-md"], ["scratch_node", "vmd"], ["scratch_node", "amber-md"],
                          ["single_step_plan", "wizard-vmd"]])
    add("graph:analysis-inputs", [["apply_template", "amber", 40, 40], ["analysis_inputs", "n4"],
                                  ["remove_node", "n3"], ["analysis_inputs", "n4"]])
    return ops


def run_graph(c, TOOLS, P):
    from mdatk.config import WorkflowParams

    wf = WorkflowParams(project_name="golden", stamp=STAMP)
    root = "/golden/plan"
    out = []
    for op in c["ops"]:
        name, args = op[0], op[1:]
        try:
            if name == "add_node":
                n = P.add_node(wf, args[0], args[1], args[2])
                res = {"id": n.id, "step": n.step}
            elif name == "connect":
                e = P.connect(wf, *args)
                res = {"edge": [e.source, e.output, e.target, e.input]}
            elif name == "remove_node":
                P.remove_node(wf, args[0])
                res = {}
            elif name == "next_step":
                res = {"next_step": P.next_step(wf)}
            elif name == "new_node_id":
                res = {"new_node_id": P.new_node_id(wf)}
            elif name == "apply_template":
                made = P.apply_template(wf, args[0], args[1], args[2])
                res = {"added": [[n.id, n.tool, n.step, n.x, n.y] for n in made]}
            elif name == "scratch_node":
                n = P.scratch_node(wf, args[0])
                res = {"id": n.id, "step": n.step, "scratch": [s.id for s in wf.scratch]}
            elif name == "single_step_plan":
                node = next(s for s in wf.scratch if s.id == args[0])
                p = P.single_step_plan(wf, node, root)
                res = {"order": p.order, "labels": [p.steps[n].label for n in p.order],
                       "issues": [[i.level, i.message] for _n, i in p.all_issues()]}
            elif name == "analysis_inputs":
                p = P.Plan(wf, root)
                top, traj = P.analysis_inputs(p, args[0])
                res = {"topology": top, "trajectories": traj}
            else:
                raise ValueError(name)
        except (KeyError, ValueError) as exc:
            res = {"error": type(exc).__name__, "message": str(exc).strip("'") if isinstance(exc, KeyError)
                   else str(exc)}
        out.append(res)
    state = {"nodes": [[n.id, n.tool, n.step, n.x, n.y] for n in wf.nodes],
             "edges": [[e.source, e.output, e.target, e.input] for e in wf.edges]}
    return {"name": c["name"], "ops": out, "final": state}


# ================================================================================================
# fixtures

L_LIB_UNITS = ("A0", "C0", "G0", "U0", "A05", "C05", "G05", "U05", "A03", "C03", "G03", "U03")
DUDEK_UNITS = ("LA", "LC", "LG", "LU", "LA5", "LC5", "LG5", "LU5", "LA3", "LC3", "LG3", "LU3")


def off_library(units):
    """A minimal OFF library text: its unit index and an atoms table per unit (what the planner reads)."""
    lines = ["!!index array str"] + [f' "{u}"' for u in units]
    for u in units:
        lines += [f"!entry.{u}.unit.atoms table  str name  str type  int typex  int resx  int flags  int seq  "
                  "int elmnt  dbl chg",
                  ' "P" "P" 0 1 131072 1 15 1.1662',
                  " \"O5'\" \"LD\" 0 1 131072 2 8 -0.4989",
                  f"!entry.{u}.unit.name single str", f' "{u}"']
    return "\n".join(lines) + "\n"


def write_fixtures(mdatk, out):
    """The fixtures: PDB models built with MDATK's builder (D, L, mixed, blank chains), an annotated model,
    secondary-structure files, residue libraries and a CP frcmod."""
    sys.path.insert(0, os.path.abspath(mdatk))
    import improvements_chem

    improvements_chem.apply_mdatk()
    from mdatk import annotate, builder

    os.makedirs(out, exist_ok=True)

    def write(name, text):
        with open(os.path.join(out, name), "w", encoding="utf-8", newline="\n") as fh:
            fh.write(text)

    d = builder.build_helix(FIXTURE_SEQUENCE, chirality="D", scheme="amber")
    lh = builder.build_helix(FIXTURE_SEQUENCE, chirality="L", scheme="amber", chains=("C", "D"), start=(7, 7))
    write("d_duplex.pdb", d)
    write("l_duplex.pdb", builder.build_helix(FIXTURE_SEQUENCE, chirality="L", scheme="amber"))
    atoms = lambda text: [ln for ln in text.splitlines() if ln.startswith(("ATOM", "HETATM", "TER"))]  # noqa: E731
    write("dl_duplex.pdb", "\n".join(atoms(d) + atoms(lh) + ["END"]) + "\n")
    # both strands with a blank chain id and numbered from 1 (as built), separated by TER (WG-06)
    blank = [line[:21] + " " + line[22:] if len(line) > 21 else line for line in atoms(d)]
    write("blank_chains.pdb", "\n".join(blank + ["END"]) + "\n")
    write("annotated.pdb", annotate.annotate(d, engine="builtin").text)
    write("hairpin.secstruct", ">hairpin\n" + TYPED_SS["dotbracket"] + "\n" + TYPED_SS["sequence"] + "\n")
    write("hairpin_no_sequence.secstruct", TYPED_SS["dotbracket"] + "\n")
    write("lrna_a0.lib", off_library(L_LIB_UNITS))
    write("lrna_dudek.lib", off_library(DUDEK_UNITS))
    write("cp_types.frcmod", "CP phosphate oxygens (test file)\nMASS\nOP 16.00\nOR 16.00\n\nNONB\n  OP 1.7493 0.2100\n"
                             "  OR 1.6837 0.1700\n\n")


# ================================================================================================
# main

def run_worker(kind, mdatk, cases_file, out_file, data_dir):
    subprocess.run([sys.executable, "-I", "-B", os.path.abspath(__file__), "--worker", kind, mdatk,
                    cases_file or "-", out_file, data_dir], check=True)
    with open(out_file, encoding="utf-8") as fh:
        return json.load(fh)


def merge(improved, raw):
    """The expected data: the improved results; where MDATK v0.0.4.1 differs, its issues and text hashes."""
    texts = dict(improved["texts"])
    by_group = {}
    for c, res, old in zip(improved["cases"], improved["results"], raw["results"]):
        entry = dict(c)
        entry.pop("group")
        entry["result"] = res
        if c["group"] != "graph":
            diff = differences(res, old, improved["texts"], raw["texts"])
            if diff:
                entry["mdatk"] = diff
        elif res != old:
            entry["mdatk"] = {"ops": old["ops"], "final": old["final"]}
        by_group.setdefault(c["group"], []).append(entry)
    return by_group, texts


def sha(text):
    return hashlib.sha256(text.encode("utf-8", "surrogateescape")).hexdigest()


D11_TEXTS = [("MDATK <VERSION>", "MDBridge <VERSION>"), ("MDATK only writes this script", "MDBridge only writes this script"),
             ("# --- settings written by MDATK ---", "# --- settings written by MDBridge ---")]


def d11(text):
    """MDATK's text with D11's program name, so that "mdatk" records the other differences only."""
    for old, new in D11_TEXTS:
        text = text.replace(old, new)
    return text


def differences(new, old, new_texts, old_texts):
    """Where MDATK v0.0.4.1 differs from the expected data beyond D11's program name: its issues, the first
    16 hex digits of the SHA-256 of each of its texts (step files as <node>/<name>, the master and the
    summary by path; null: a file MDATK does not write), the steps it skips; where it crashes, its exception."""
    if "crash" in old:
        return {"crash": old["crash"]}
    out = {}
    old_issues = [[n, level, d11(m), keys] for n, level, m, keys in old["issues"]]
    if new["issues"] != old_issues:
        out["issues"] = old_issues

    def files(res, store):
        got = {}
        for nid, step in res["steps"].items():
            for f in step.get("files", []):
                got[f"{nid}/{f['name']}"] = d11(store[f["text"]])
        w = res["write"]
        for rel in (w["master"], w["summary"]):
            if rel:
                got[rel] = d11(store[w["written"][rel]["text"]])
        return got

    a, b = files(new, new_texts), files(old, old_texts)
    texts = {}
    for key in sorted(set(a) | set(b)):
        if a.get(key) != b.get(key):
            texts[key] = sha(b[key])[:16] if key in b else None
    if texts:
        out["texts"] = texts
    if new["write"]["skipped"] != {k: d11(v) for k, v in old["write"]["skipped"].items()}:
        out["skipped"] = old["write"]["skipped"]
    return out


def main(argv):
    if argv[:1] == ["--worker"]:
        kind, mdatk, cases_file, out_file, data_dir = argv[1:6]
        worker(kind, mdatk, None if cases_file == "-" else cases_file, out_file, data_dir)
        return 0
    if argv[:1] == ["--fixtures"]:
        write_fixtures(argv[1], argv[2])
        return 0
    if len(argv) != 2:
        print(__doc__.split("\n\n")[1], file=sys.stderr)
        return 2
    mdatk, out = os.path.abspath(argv[0]), os.path.abspath(argv[1])
    fixtures = os.path.join(out, "fixtures")
    for sub in ("fixtures", "cases", "texts"):
        shutil.rmtree(os.path.join(out, sub), ignore_errors=True)
    subprocess.run([sys.executable, "-I", "-B", os.path.abspath(__file__), "--fixtures", mdatk, fixtures], check=True)
    # demo_dry.parm7 (topology fixtures of tests/data) is read in place, as ../topology/demo_dry.parm7
    tmp = tempfile.mkdtemp(prefix="workflow-cases-")
    try:
        data_dir = os.path.dirname(out)
        improved = run_worker("improved", mdatk, None, os.path.join(tmp, "improved.json"), data_dir)
        cases_file = os.path.join(tmp, "cases.json")
        with open(cases_file, "w", encoding="utf-8") as fh:
            json.dump(improved["cases"], fh)
        raw = run_worker("mdatk", mdatk, cases_file, os.path.join(tmp, "raw.json"), data_dir)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    by_group, texts = merge(improved, raw)
    os.makedirs(os.path.join(out, "cases"))
    os.makedirs(os.path.join(out, "texts"))
    used = set()
    for group, cases in sorted(by_group.items()):
        for c in cases:
            r = c["result"]
            for step in r.get("steps", {}).values():
                used.update(f["text"] for f in step.get("files", []))
            used.update(f["text"] for f in r.get("write", {}).get("written", {}).values())
        with open(os.path.join(out, "cases", f"{group}.json"), "w", encoding="utf-8", newline="\n") as fh:
            fh.write('{"stamp": "%s", "root": "%s", "texts": "../texts/<id>.txt", "cases": [\n' % (STAMP, ROOT))
            fh.write(",\n".join(json.dumps(c, ensure_ascii=False, separators=(",", ":")) for c in cases))
            fh.write("\n]}\n")
    for key in sorted(used):
        with open(os.path.join(out, "texts", key + ".txt"), "w", encoding="utf-8", errors="surrogateescape",
                  newline="") as fh:
            fh.write(texts[key])
    n = sum(len(c) for c in by_group.values())
    print(f"{n} cases in {len(by_group)} groups, {len(used)} texts -> {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
