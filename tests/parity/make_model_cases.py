"""Write the test vectors of MDBridge's parameter model (core/project.hpp, mdatk/config.py).

    python3 -I -B tests/parity/make_model_cases.py MDATK_CHECKOUT OUT_DIR

Run with Python 3.13 on Linux (repr() and the Unicode data depend on the version), numpy and scipy
importable (density_problems imports mdatk.density).
-I / -B: the MDATK checkout is put on sys.path explicitly and gets no byte code. Writes into OUT_DIR
(tests/data/project-model):

* ``from_dict_cases.json``  Project.from_dict(input, base_dir): repr() of to_dict() (or of the entry at
                            "result_path"), or the error message; the warnings logged on the way.
                            "path" + "value" cases set one entry (a list index 0 stands for a one-item
                            list of objects), "input" cases give the whole data.
* ``validate_cases.json``   validate(), density_problems(), density_active(), layout_warnings() and
                            group_names() of Project.from_dict(input), run in a folder holding the
                            files of "files" (empty).
* ``config_doc.json``       projectfile.config_doc of every path of the model (list indices left out).
* ``default_project.json``  repr() of Project().to_dict() and its keys in order.

Tagged values ("value" / "input"): {"t": "none" | "bool" | "int" | "float" | "str" | "datetime" | "date" |
"time" | "array" | "table", ...}, the format of make_writer_vectors.py (tests/core/test_support.hpp).
"""

from __future__ import annotations

import dataclasses
import datetime as dt
import json
import logging
import math
import os
import sys
import tempfile
import typing

if len(sys.argv) != 3:
    sys.exit("usage: python3 -I -B tests/parity/make_model_cases.py MDATK_CHECKOUT OUT_DIR")
MDATK, OUT = (os.path.abspath(a) for a in sys.argv[1:3])
sys.dont_write_bytecode = True
sys.path.insert(0, MDATK)

from mdatk import config as C  # noqa: E402
from mdatk import projectfile as PF  # noqa: E402

os.makedirs(OUT, exist_ok=True)

NAN, INF = float("nan"), float("inf")


def dump(name, data):
    """A JSON object whose lists hold one item per line."""
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
    if isinstance(v, dt.date):
        return {"t": "date", "iso": v.isoformat()}
    if isinstance(v, dt.time):
        return {"t": "time", "iso": v.isoformat()}
    if isinstance(v, (list, tuple)):
        return {"t": "array", "v": [enc(x) for x in v]}
    if isinstance(v, dict):
        return {"t": "table", "v": [[k, enc(x)] for k, x in v.items()]}
    raise TypeError(type(v))


class Capture(logging.Handler):
    def __init__(self):
        super().__init__(logging.DEBUG)
        self.messages = []

    def emit(self, record):
        self.messages.append(record.getMessage())


CAPTURE = Capture()
logging.getLogger("mdatk").addHandler(CAPTURE)
logging.getLogger("mdatk").setLevel(logging.DEBUG)
logging.getLogger("mdatk").propagate = False


def from_dict(data, base_dir=None):
    """(project or None, error or None, warnings)"""
    CAPTURE.messages = []
    try:
        return C.Project.from_dict(data, base_dir=base_dir), None, list(CAPTURE.messages)
    except (ValueError, OverflowError) as exc:
        return None, str(exc), list(CAPTURE.messages)


def at(plain, path):
    for key in path:
        plain = plain[key]
    return plain


# ------------------------------------------------------------------------------ the model's paths


def unwrap(tp):
    while True:
        origin, args = typing.get_origin(tp), typing.get_args(tp)
        if origin is typing.Union:
            tp = [a for a in args if a is not type(None)][0]
        elif origin in (list, typing.List) and args:
            tp = args[0]
        else:
            return tp


def model_paths(cls=C.Project, prefix=()):
    """(path, type) of every field, depth first in declaration order."""
    hints = typing.get_type_hints(cls)
    for f in dataclasses.fields(cls):
        path = prefix + (f.name,)
        yield path, hints[f.name]
        sub = unwrap(hints[f.name])
        if dataclasses.is_dataclass(sub):
            yield from model_paths(sub, path)


def leaf_paths(cls=C.Project, prefix=()):
    """Input paths of every non-object field; a list of objects contributes index 0."""
    hints = typing.get_type_hints(cls)
    for f in dataclasses.fields(cls):
        tp = hints[f.name]
        sub = unwrap(tp)
        if dataclasses.is_dataclass(sub):
            list_of = typing.get_origin(tp) in (list, typing.List)
            yield from leaf_paths(sub, prefix + ((f.name, 0) if list_of else (f.name,)))
        else:
            yield prefix + (f.name,), tp


def nest(path, value):
    """The input that sets `path` (an int: a one-item list of objects)."""
    for key in reversed(path):
        value = [value] if isinstance(key, int) else {key: value}
    return value


# ------------------------------------------------------------------------------------- from_dict

# a value of every TOML type, and the edge cases of the numeric checks
SAMPLES = [
    None, True, False, 0, 7, -3, 2**63 - 1, -2**63, 2.5, 4.0, -0.0, 1e-07, 1e16, NAN, INF, -INF,
    "", "text", "4.5", "é \"q\" \n", dt.datetime(2026, 10, 15, 14, 32, 5),
    dt.datetime(2026, 10, 15, 14, 32, 5, 250000, tzinfo=dt.timezone(dt.timedelta(hours=2))),
    dt.date(2026, 10, 15), dt.time(14, 32, 5), [], [1, "a", 2.5], [[1], {"k": None}],
    {}, {"k": 1, "n": None},
]
# floats an int field takes in Python but MDBridge rejects (big integers; see project.hpp)
BEYOND_INT64 = (2.0**63, -2.0**64, 1e300)


def leaf_cases():
    out = []
    for path, _tp in leaf_paths():
        for v in SAMPLES:
            p, err, warns = from_dict(nest(path, v))
            case = {"path": list(path), "value": enc(v)}
            if err is not None:
                case["error"] = err
            else:
                case["result"] = repr(at(p.to_dict(), path))
            if warns:
                case["warnings"] = warns
            out.append(case)
    return out


def float_int_cases():
    out = []
    for path, tp in leaf_paths():
        if unwrap(tp) is not int:
            continue
        for v in (3.0, -0.0, 2.0**62, -2.0**63, 0.5, -1.5, 9007199254740993.0):
            p, err, warns = from_dict(nest(path, v))
            case = {"path": list(path), "value": enc(v)}
            if err is not None:
                case["error"] = err
            else:
                case["result"] = repr(at(p.to_dict(), path))
            out.append(case)
    return out


def whole_cases():
    """Hand-made inputs: (name, data, base_dir, result_path)."""
    big_options = {"a": None, "b": {"c": None, "d": [None, {"e": None, "f": 1}], "g": {}},
                   "mode": "frequency", "when": dt.datetime(2026, 1, 2, 3, 4, 5)}
    node = {"id": "n1", "tool": "input-file", "step": 0,
            "settings": {"path": "maps/m.mrc", "kind": "map", "x": None, "deep": {"y": None, "z": 1}},
            "files": {"structure": " a.pdb ; ;b.pdb ", "none": None, "num": 5}}
    base_node = dict(node, files={"structure": " a.pdb ; ;/b.pdb ", "blank": " ; ", "empty": ""})
    return [
        ("empty", {}, None, ()),
        ("meta keys are dropped", {"format": "mdatk-project", "format_version": 1, "mdatk_version": "0.0.4.1",
                                   "saved": dt.datetime(2026, 10, 15, 14, 32), "log": [{"step": "x"}],
                                   "log_file": "p.log.toml", "notes": "kept"}, None, ()),
        ("meta keys only at the top", {"inputs": {"format": 1, "log": 2}}, None, ("inputs",)),
        ("unknown keys everywhere", {"bogus": 1, "inputs": {"nope": 2, "topology": "t"},
                                     "analyses": {"contacts": {"far": True}},
                                     "groups": [{"name": "A", "residues": "1", "extra": 1}],
                                     "layout": {"panels": [{"type": "rmsd", "what": 1}]},
                                     "workflow": {"nodes": [{"id": "x", "port": 1}],
                                                  "edges": [{"source": "x", "weight": 1}]},
                                     "é\n'\"": 1, "": 2}, None, ()),
        ("alias input", {"input": {"topology": "t.parm7"}}, None, ("inputs",)),
        ("alias trajectory string", {"inputs": {"trajectory": "t.nc"}}, None, ("inputs",)),
        ("alias trajectory list", {"inputs": {"trajectory": ["a.nc", "b.nc"]}}, None, ("inputs",)),
        ("trajectories string", {"inputs": {"trajectories": "a.nc"}}, None, ("inputs",)),
        ("trajectories empty string", {"inputs": {"trajectories": ""}}, None, ("inputs",)),
        ("trajectories None", {"inputs": {"trajectories": None}}, None, ("inputs",)),
        ("trajectory then trajectories", {"inputs": {"trajectory": "a.nc", "trajectories": ["b.nc"]}},
         None, ("inputs",)),
        ("trajectories then trajectory", {"inputs": {"trajectories": ["b.nc"], "trajectory": "a.nc"}},
         None, ("inputs",)),
        ("input then inputs, no merge", {"input": {"topology": "a", "step": 3}, "inputs": {"start": 2}},
         None, ("inputs",)),
        ("inputs then input, no merge", {"inputs": {"start": 2}, "input": {"topology": "a"}}, None, ("inputs",)),
        ("alias inside the wrong table", {"inputs": {"input": 1}, "trajectory": "x"}, None, ()),
        ("bad value under an alias", {"input": {"trajectory": 5}}, None, ()),
        ("bad later duplicate", {"input": {"topology": "a"}, "inputs": 3}, None, ()),
        ("bad earlier duplicate", {"input": 3, "inputs": {"topology": "a"}}, None, ()),
        ("groups replace the defaults", {"groups": [{"name": "X"}]}, None, ("groups",)),
        ("no groups", {"groups": []}, None, ("groups",)),
        ("groups None", {"groups": None}, None, ("groups",)),
        ("rmsd targets replace the defaults", {"analyses": {"rmsd": {"targets": [{"label": "L"}]}}},
         None, ("analyses", "rmsd")),
        ("panels replace the defaults", {"layout": {"panels": [{"type": "timeseries"}]}}, None, ("layout",)),
        ("float field given an int", {"analyses": {"contacts": {"cutoff": 3}}}, None,
         ("analyses",)),
        ("object item not a table", {"groups": [{"name": "A"}, 5]}, None, ()),
        ("object item a list", {"layout": {"panels": [[1]]}}, None, ()),
        ("object None", {"analyses": None}, None, ()),
        ("object a string", {"output": {"video": "v.mp4"}}, None, ()),
        ("list a string", {"visual": {"rotation": "0 0 0"}}, None, ()),
        ("list a table", {"output": {"charts": {"formats": {"png": True}}}}, None, ()),
        ("list items", {"visual": {"rotation": [1, 2.5, -0.0]}}, None, ("visual", "rotation")),
        ("list item wrong", {"visual": {"rotation": [1, "x", 3]}}, None, ()),
        ("list item bool", {"visual": {"rotation": [1, 2, True]}}, None, ()),
        ("anchor", {"density": {"anchor": [[1, 2, 3], [4.5, 5, 6], [], None]}}, None, ("density", "anchor")),
        ("anchor item wrong", {"density": {"anchor": [[1, 2, 3], [4, "x"]]}}, None, ()),
        ("anchor item a number", {"density": {"anchor": [[1, 2, 3], 4]}}, None, ()),
        ("str list items coerced", {"output": {"charts": {"series": [1, 2.5, True, None, [1, "a"], {"k": "v"},
                                                                      dt.date(2026, 1, 2)]}}},
         None, ("output", "charts", "series")),
        ("options keep values, drop None at every depth",
         {"layout": {"panels": [{"type": "contact_map", "options": big_options}]}}, None, ("layout", "panels")),
        ("options not a table", {"layout": {"panels": [{"options": [1]}]}}, None, ()),
        ("options None", {"layout": {"panels": [{"options": None}]}}, None, ()),
        ("files and settings unchecked", {"workflow": {"nodes": [node], "scratch": [node]}}, None, ("workflow",)),
        ("files not a table", {"workflow": {"nodes": [{"files": "a.pdb"}]}}, None, ()),
        ("optional None", {"inputs": {"stop": None, "dt_ps": None}, "analyses": {"energy": {"target_value": None}}},
         None, ()),
        ("optional set", {"inputs": {"stop": 10, "dt_ps": 2}, "cleaning": {"stop": 5.0},
                          "analyses": {"energy": {"target_value": -3}}}, None, ()),
        ("first error wins, in input order", {"output": {"video": {"fps": "x"}}, "analyses": {"contacts": {
            "cutoff": "far"}}}, None, ()),
        ("first error wins, reversed", {"analyses": {"contacts": {"cutoff": "far"}}, "output": {"video": {
            "fps": "x"}}}, None, ()),
        ("warnings before the error", {"bogus": 1, "inputs": {"what": 1, "step": "x"}, "late": 2}, None, ()),
        ("error in a later list item", {"groups": [{"name": "A"}, {"name": "B", "color": None},
                                                   {"name": 3, "residues": [1]}, {"bogus": 1}]}, None, ()),
        ("enabled given 1", {"analyses": {"rmsd": {"enabled": 1}}}, None, ()),
        ("cutoff given far", {"analyses": {"contacts": {"cutoff": "far"}}}, None, ()),
        ("density round trip", {"density": {"map_path": "m.mrc", "fit_frame": "first", "transform": [1.0, 2.0]}},
         None, ("density",)),
        ("full project", rich_project().to_dict(), None, ()),
        # base_dir: relative paths resolved
        ("base dir", {"inputs": {"topology": "top.parm7", "trajectory": "../t.nc",
                                 "external": [{"path": "data/x.dat"}, {"path": "/abs/y.dat"}]},
                      "cleaning": {"cpptraj_path": "bin/cpptraj", "reference_file": "r.pdb", "output_dir": "c"},
                      "output": {"out_dir": "out", "cache_dir": "./cache/"},
                      "density": {"map_path": "m.mrc"}, "workflow": {"out_dir": "wf", "nodes": [base_node]}},
         "/base/proj", ()),
        ("base dir, defaults", {}, "/base/proj", ()),
        ("base dir relative", {"inputs": {"topology": "t"}, "output": {"out_dir": "../o"}}, "rel/dir", ()),
        ("base dir empty: no resolution", {"inputs": {"topology": "t"}}, "", ()),
        ("base dir, cpptraj on PATH", {"cleaning": {"cpptraj_path": "cpptraj"}}, "/b", ("cleaning",)),
        ("base dir, absolute and dotted paths", {"inputs": {"topology": "/x/../t", "trajectories": ["a/./b", "."]},
                                                 "output": {"out_dir": ""}}, "/b/c/", ()),
        ("base dir, unknown step type and string settings",
         {"workflow": {"nodes": [{"tool": "amber-build", "settings": {"path": "p.pdb"}, "files": {"a": "b"}},
                                 {"tool": "input-file", "settings": {"path": 3}},
                                 {"tool": "input-file", "settings": {"path": ""}},
                                 {"tool": "input-file", "settings": {"kind": "pdb"}}]}}, "/b", ("workflow",)),
    ]


def rich_project():
    p = C.Project()
    p.notes = "Aim: compare\nthree \"quoted\" runs"
    p.inputs.topology = "/data/sys.parm7"
    p.inputs.trajectories = ["/data/a.nc", "/data/b.nc"]
    p.inputs.stop = 500
    p.inputs.dt_ps = 2.5
    p.inputs.external = [C.ExternalSeries("/data/e.dat", "E", "0", "ns", ["1", "rmsd"], ",")]
    p.cleaning.enabled = True
    p.cleaning.stop = 100
    p.analyses.contacts.cutoff = 5.25
    p.analyses.energy.target_value = -42.5
    p.visual.representations = [C.RepSpec("licorice", ["Aptamer"], "element", 0.5, 0.8, True)]
    p.visual.rotation = [10.0, -20.0, 30.5]
    p.density.map_path = "/data/m.mrc"
    p.density.transform = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.5, -2.0, 3.0]
    p.density.anchor = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]
    p.layout.panels.append(C.PanelSpec("contact_map", ["contacts/map"], "Map", {"mode": "frequency"}))
    p.output.charts.formats = ["png", "svg"]
    p.workflow.nodes = [C.PlanNode("n1", "input-file", 0, "In", 10.0, 20.0, {"path": "/m.pdb", "kind": "pdb"},
                                   {"structure": "/m.pdb"})]
    p.workflow.edges = [C.PlanEdge("n1", "file", "n2", "structure")]
    return p


def from_dict_cases():
    cases = []
    for name, data, base, result_path in whole_cases():
        p, err, warns = from_dict(data, base)
        case = {"name": name, "input": enc(data), "base_dir": base, "result_path": list(result_path)}
        if err is not None:
            case["error"] = err
        else:
            case["result"] = repr(at(p.to_dict(), result_path))
        case["warnings"] = warns
        cases.append(case)
    return {"how": "Project.from_dict(input, base_dir); repr(to_dict()) at result_path, or str(error); "
                   "warnings logged by mdatk.config",
            "cases": cases, "leaf_cases": leaf_cases(), "float_int_cases": float_int_cases(),
            "beyond_int64": [enc(v) for v in BEYOND_INT64]}


# -------------------------------------------------------------------------------------- validate

FILES = ["top.parm7", "a.nc", "b.nc", "map.mrc"]
ROT90 = [0.0, -1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0]
IDENT = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0]
DENS = {"map_path": "map.mrc"}


def dens(**kw):
    return {"density": dict(DENS, **kw)}


def validate_inputs():
    reps = [{"style": "cartoon"}, {"style": "spheres", "color": "rainbow", "groups": ["Aptamer", "Nope", ""]},
            {"size": 0.0}, {"opacity": 0.0}, {"opacity": 1.0}, {"opacity": 1.5}, {"size": NAN}, {"opacity": NAN},
            {"size": -1.0, "opacity": -1.0}, {"style": "trace", "color": "nucleotide", "size": 1e-9}]
    out = [
        ("default", {}),
        ("inputs ok", {"inputs": {"topology": "top.parm7", "trajectories": ["a.nc", "b.nc"]}}),
        ("inputs missing", {"inputs": {"topology": "nope.parm7", "trajectories": ["a.nc", "x.nc", "y.nc", ""]}}),
        ("topology is a folder", {"inputs": {"topology": ".", "trajectories": ["sub"]}}),
        ("step 0", {"inputs": {"step": 0}}),
        ("step -1", {"inputs": {"step": -1}}),
        ("step 1", {"inputs": {"step": 1}}),
        ("duplicate groups", {"groups": [{"name": "A", "residues": "1"}, {"name": "A", "residues": "2"}]}),
        ("group names", {"groups": [{"name": "", "residues": "1"}, {"name": "  ", "residues": "1"},
                                    {"name": "Bad name", "residues": "1-3"}, {"name": "ok.name-1_", "selection": "x"},
                                    {"name": "café", "residues": "1"}, {"name": "a\n", "residues": "1"},
                                    {"name": "x/y", "residues": " "}, {"name": "E", "residues": "\t",
                                                                       "selection": "　"},
                                    {"name": " ", "residues": "1"}, {"name": "'q\"", "residues": "1"}]}),
        ("no groups", {"groups": []}),
        ("unknown analysis groups", {"analyses": {
            "contacts": {"group_a": "X", "group_b": ""}, "hbonds": {"group_a": "Aptamer", "group_b": "Y"},
            "mg": {"group_a": "M1", "group_b": "M2"}, "stacking": {"group_a": "S", "group_b": "aptamer"},
            "energy": {"group_a": "E", "group_b": "Target"}}}),
        ("disabled analyses are not checked", {"analyses": {
            "contacts": {"enabled": False, "group_a": "X"}, "hbonds": {"enabled": False, "group_b": "Y"},
            "mg": {"enabled": False, "group_a": "M1"}, "stacking": {"enabled": False, "group_a": "S"},
            "energy": {"enabled": False, "group_a": "E"}}}),
        ("contacts cutoff 0", {"analyses": {"contacts": {"cutoff": 0}}}),
        ("contacts cutoff -1", {"analyses": {"contacts": {"cutoff": -1.0}}}),
        ("contacts cutoff nan", {"analyses": {"contacts": {"cutoff": NAN}}}),
        ("contacts cutoff tiny", {"analyses": {"contacts": {"cutoff": 5e-324}}}),
        ("mg inner > outer", {"analyses": {"mg": {"inner_cutoff": 5.0, "outer_cutoff": 4.5}}}),
        ("mg inner = outer", {"analyses": {"mg": {"inner_cutoff": 4.5}}}),
        ("mg nan", {"analyses": {"mg": {"inner_cutoff": NAN}}}),
        ("mg inf", {"analyses": {"mg": {"inner_cutoff": INF, "outer_cutoff": INF}}}),
        ("modes", {"visual": {"color_mode": "rainbow", "glow_mode": "sparkle"}}),
        ("modes ok", {"visual": {"color_mode": "chain", "glow_mode": "none"}}),
        ("modes case", {"visual": {"color_mode": "Groups", "glow_mode": "FADE"}}),
        ("representations", {"visual": {"representations": reps}}),
        ("renderer", {"visual": {"renderer": "opengl"}}),
        ("renderers ok", {"visual": {"renderer": "vtk"}}),
        ("zoom low", {"visual": {"zoom": 0.04}}),
        ("zoom bounds", {"visual": {"zoom": 0.05}}),
        ("zoom high", {"visual": {"zoom": 50.0001}}),
        ("zoom 50", {"visual": {"zoom": 50}}),
        ("zoom nan", {"visual": {"zoom": NAN}}),
        ("base pairs", {"analyses": {"basepairs": {"engine": "x3dna", "stride": 0, "cutoff": 1.9}}}),
        ("base pairs high", {"analyses": {"basepairs": {"engine": "barnaba", "stride": -5, "cutoff": 5.1}}}),
        ("base pairs bounds", {"analyses": {"basepairs": {"stride": 1, "cutoff": 2}}}),
        ("base pairs nan", {"analyses": {"basepairs": {"cutoff": NAN}}}),
        ("rep rest and backbone", {"visual": {"rep_rest": "none", "backbone_atom": "CA"}}),
        ("backbone C4'", {"visual": {"backbone_atom": "C4'", "rep_rest": "hide"}}),
        ("proximity", {"visual": {"proximity_cutoff": 8.5}}),
        ("proximity equal", {"visual": {"proximity_cutoff": 8}}),
        ("proximity nan", {"visual": {"proximity_cutoff": NAN}}),
        ("proximity radius nan", {"analyses": {"proximity_search_radius": NAN}}),
        ("fps 0", {"output": {"video": {"fps": 0}}}),
        ("fps -1", {"output": {"video": {"fps": -1}}}),
        ("density ok", dens()),
        ("density missing map", {"density": {"map_path": "gone.mrc"}}),
        ("density missing map, hidden", {"density": {"map_path": "gone.mrc", "show": False}}),
        ("density missing map, no structure panel", {"density": {"map_path": "gone.mrc"},
                                                     "layout": {"panels": [{"type": "timeseries"}]}}),
        ("density missing map, no panels", {"density": {"map_path": "gone.mrc"}, "layout": {"panels": []}}),
        ("density no map, everything wrong", {"density": {"resolution": -1, "opacity": 2, "angular_step": 1,
                                                          "search": "x", "transform": [1], "fit_frame": "q"}}),
        ("density resolution 0", dens(resolution=0)),
        ("density resolution nan", dens(resolution=NAN)),
        ("density opacity", dens(opacity=-0.1)),
        ("density opacity high", dens(opacity=1.1)),
        ("density opacity bounds", dens(opacity=0)),
        ("density opacity 1", dens(opacity=1)),
        ("density opacity nan", dens(opacity=NAN)),
        ("density angular step low", dens(angular_step=4.99)),
        ("density angular step high", dens(angular_step=46)),
        ("density angular step bounds", dens(angular_step=5)),
        ("density angular step 45", dens(angular_step=45)),
        ("density angular step nan", dens(angular_step=NAN)),
        ("density search", dens(search="Global")),
        ("density search manual", dens(search="manual")),
        ("transform identity", dens(transform=IDENT + [1, 2, 3])),
        ("transform rotation", dens(transform=ROT90 + [0, 0, 0])),
        ("transform short", dens(transform=[1.0, 2.0])),
        ("transform long", dens(transform=IDENT + [0, 0, 0, 0])),
        ("transform inf", dens(transform=[1, 0, 0, 0, 1, 0, 0, 0, 1, INF, 0, 0])),
        ("transform nan", dens(transform=[NAN] + IDENT[1:] + [0, 0, 0])),
        ("transform nan and short", dens(transform=[NAN])),
        ("transform scaled", dens(transform=[2 * x for x in IDENT] + [0, 0, 0])),
        ("transform near rotation", dens(transform=[1.005, 0, 0, 0, 0.995, 0, 0, 0, 1, 0, 0, 0])),
        ("transform off rotation", dens(transform=[1.02, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0])),
        ("transform reflection", dens(transform=[-1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0])),
        ("transform reflection rotated", dens(transform=[0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0])),
        ("transform singular", dens(transform=[1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0])),
        ("transform zeros", dens(transform=[0.0] * 12)),
        ("transform sheared", dens(transform=[1, 0.005, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0])),
        ("transform sheared more", dens(transform=[1, 0.03, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0])),
        ("transform rounded rotation", dens(transform=[0.707107, -0.707107, 0, 0.707107, 0.707107, 0,
                                                       0, 0, 1, 5, 5, 5])),
        ("transform big translation", dens(transform=IDENT + [1e300, -1e300, 0])),
    ]
    for frame in ("first", "", " ", "ref", " REF ", "Ref", "12", " 7 ", "²", "١٢", "-1", "1.5",
                  "+3", "1 2", "0", "①", "½", "rеf", " ref\n"):
        out.append((f"fit frame {frame!r}", dens(fit_frame=frame)))
    out += [
        ("two structure panels", {"layout": {"panels": [{"type": "structure"}, {"type": "structure"}]}}),
        ("three structure panels", {"layout": {"panels": [{}, {}, {"type": "Structure"}, {"type": "structure"}]}}),
        ("everything", {
            "inputs": {"topology": "nope", "trajectories": ["x.nc"], "step": 0},
            "groups": [{"name": "A b", "residues": ""}, {"name": "A b"}],
            "analyses": {"contacts": {"cutoff": 0, "group_a": "Q"}, "mg": {"inner_cutoff": 9},
                         "basepairs": {"engine": "e", "stride": 0, "cutoff": 9}},
            "visual": {"color_mode": "c", "glow_mode": "g", "representations": reps, "renderer": "r",
                       "zoom": 100, "rep_rest": "r", "backbone_atom": "b", "proximity_cutoff": 99},
            "density": {"map_path": "gone.mrc", "resolution": 0, "opacity": 9, "angular_step": 99,
                        "search": "s", "transform": [1], "fit_frame": "f"},
            "layout": {"panels": [{}, {}]},
            "output": {"video": {"fps": 0}}}),
    ]
    return out


def validate_cases():
    old = os.getcwd()
    with tempfile.TemporaryDirectory() as tmp:
        os.chdir(tmp)
        for name in FILES:
            open(name, "w").close()
        os.mkdir("sub")
        cases = []
        for name, data in validate_inputs():
            p, err, _warns = from_dict(data)
            assert err is None, (name, err)
            cases.append({"name": name, "input": enc(data), "problems": p.validate(),
                          "density_problems": p.density_problems(), "density_active": p.density_active(),
                          "layout_warnings": C.layout_warnings(p), "group_names": p.group_names()})
        os.chdir(old)
    return {"how": "Project.from_dict(input) in a folder holding the files of 'files' (empty) and the folder "
                   "'sub': validate(), density_problems(), density_active(), layout_warnings(), group_names()",
            "files": FILES, "folders": ["sub"], "cases": cases}


# ------------------------------------------------------------------------------------ config_doc


def config_doc_cases():
    paths = []
    for path, _tp in model_paths():
        paths.append({"path": ".".join(path), "doc": PF.config_doc(path)})
    return {"how": "projectfile.config_doc(path) of every path of the model (list indices left out), "
                   "in declaration order", "paths": paths}


def default_project():
    plain = C.Project().to_dict()

    def keys(obj, prefix=""):
        out = []
        if isinstance(obj, dict):
            for k, v in obj.items():
                out.append(prefix + k)
                out += keys(v, prefix + k + ".")
        elif isinstance(obj, list):
            for i, v in enumerate(obj):
                out += keys(v, f"{prefix[:-1]}[{i}].")
        return out

    return {"how": "repr(Project().to_dict()); its keys in order", "repr": repr(plain), "keys": keys(plain)}


def main() -> int:
    dump("from_dict_cases.json", from_dict_cases())
    dump("validate_cases.json", validate_cases())
    dump("config_doc.json", config_doc_cases())
    dump("default_project.json", default_project())
    return 0


if __name__ == "__main__":
    sys.exit(main())
