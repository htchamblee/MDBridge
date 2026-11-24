"""What MDATK's citation functions (mdatk/citations.py) give, for tests/core/test_citations.cpp.

    python3 -I -B make_citation_cases.py MDATK_CHECKOUT OUT_DIR

Imports MDATK (a checkout at v0.0.4.1; nothing is written into it) and writes
OUT_DIR/citation_cases.json:

  references, resources   REFERENCES and RESOURCES in their order, every field
  unique, keys_for,       unique / keys_for / references (the keys) / bibtex / text_lines /
  references, bibtex,     grouped of many id and key lists (single, pairs, all, duplicates,
  text_lines, grouped     unknown and empty entries, random subsets)
  comment_lines           comment_lines(ids, prefix, width, heading): the lines or textwrap's
                          ValueError, for many prefixes, widths and headings
  pipeline                pipeline_resources of projects made with Project.from_dict(plain)
  project                 project_resources(project, plan) with stand-in plans: "steps" lists
                          {generates, resources} in plan order; "generating" the resource lists
                          of the steps whose tool generates a script (what MDBridge's caller
                          passes)

The output is deterministic (fixed random seed).
"""

from __future__ import annotations

import json
import os
import random
import sys

if len(sys.argv) != 3:
    sys.exit("usage: make_citation_cases.py MDATK_CHECKOUT OUT_DIR")
sys.path.insert(0, os.path.abspath(sys.argv[1]))

from mdatk import citations as C  # noqa: E402
from mdatk import config  # noqa: E402

RNG = random.Random(20261124)
IDS = list(C.RESOURCES)
KEYS = list(C.REFERENCES)


def id_lists():
    out = [[], [""], ["nope"], IDS, list(reversed(IDS)), IDS + IDS, ["numpy", "numpy", "", "scipy"],
           ["coot", "phenix_elbow", "phenix", "acedrg"], ["case2023", "ambertools"]]
    out += [[i] for i in IDS]
    out += [[a, b] for a, b in zip(IDS, IDS[1:] + IDS[:1])]
    for _ in range(40):
        pool = IDS + ["", "unknown", "mdanalysis"]
        out.append([RNG.choice(pool) for _ in range(RNG.randint(1, 12))])
    return out


def key_lists():
    out = [[], [""], ["nope"], KEYS, list(reversed(KEYS)), KEYS + KEYS, ["harris2020", "", "harris2020"]]
    out += [[k] for k in KEYS]
    for _ in range(30):
        pool = KEYS + ["", "unknown", "MDAnalysis"]
        out.append([RNG.choice(pool) for _ in range(RNG.randint(1, 10))])
    return out


def comment_cases():
    out = []
    lists = [[], ["nope"], ["mdanalysis"], ["ffmpeg"], ["ol3", "water_opc", "mg_opc"],
             ["coot", "coot", "acedrg"], IDS, ["numpy", "unknown", "scipy", "numpy"],
             ["tachyon", "vmd"], ["pdb-4jrd", "arnott-arna"]]
    for _ in range(12):
        lists.append(RNG.sample(IDS, RNG.randint(1, 8)))
    prefixes = ["# ", "", "// ", "% ", "#  ", "§ ", ";; ", "REMARK   "]
    widths = [100, 80, 60, 40, 20, 12, 7, 6, 5, 3, 2, 1, 0, -4, 200]
    headings = [None, "", "Please cite", "Cite", "References for this step, which are many words long"]
    for ids in lists:
        out.append({"ids": ids, "prefix": "# ", "width": 100, "heading": "Please cite"})
    for ids in lists[:6] + lists[-4:]:
        for prefix in prefixes:
            for width in widths:
                out.append({"ids": ids, "prefix": prefix, "width": width,
                            "heading": RNG.choice(headings)})
    for ids in lists[:4]:
        for heading in headings:
            out.append({"ids": ids, "prefix": "# ", "width": 70, "heading": heading})
    for c in out:
        try:
            c["lines"] = C.comment_lines(c["ids"], prefix=c["prefix"], width=c["width"],
                                         heading=c["heading"])
        except ValueError as e:
            c["error"] = str(e)
    return out


def panels(*specs):
    return [{"type": t, "options": o} if o is not None else {"type": t} for t, o in specs]


def project_variants():
    out = [{}]
    for cleaning in (False, True):
        for video in (False, True):
            for renderer in ("auto", "software", "vtk", "VTK"):
                out.append({"cleaning": {"enabled": cleaning}, "output": {"video": {"enabled": video}},
                            "visual": {"renderer": renderer}})
    panel_sets = [
        [],
        panels(("structure", None)),
        panels(("secondary", None)),
        panels(("secondary", {})),
        panels(("secondary", {"source": "dotbracket"})),
        panels(("secondary", {"source": "dynamic"})),
        panels(("secondary", {"source": "static"})),
        panels(("secondary", {"source": "Dotbracket"})),
        panels(("secondary", {"source": 3})),
        panels(("secondary", {"source": True})),
        panels(("secondary", {"source": ["dotbracket"]})),
        panels(("secondary", {"mode": "dotbracket"})),
        panels(("secondary", {"source": "dotbracket"}), ("secondary", {"source": "dynamic"})),
        panels(("secondary", {"source": "dotbracket"}), ("structure", None)),
        panels(("timeseries", {"source": "dynamic"}), ("contact_map", None)),
        panels(("Secondary", None)),
    ]
    for engine in ("builtin", "barnaba", "Barnaba", "dssr"):
        for ps in panel_sets:
            out.append({"analyses": {"basepairs": {"engine": engine}}, "layout": {"panels": ps}})
    out.append({"cleaning": {"enabled": True}, "output": {"video": {"enabled": True}},
                "visual": {"renderer": "vtk"}, "analyses": {"basepairs": {"engine": "barnaba"}},
                "layout": {"panels": panels(("secondary", None))}})
    out.append({"cleaning": {"enabled": True}, "output": {"video": {"enabled": False}},
                "visual": {"renderer": "vtk"}, "analyses": {"basepairs": {"engine": "barnaba"}},
                "layout": {"panels": panels(("structure", None), ("secondary", {"source": "x"}))}})
    return out


class Tool:
    def __init__(self, generates, resources):
        self.generates = generates
        self._resources = resources

    def resources(self, settings):
        return list(self._resources)


class Step:
    def __init__(self, tool):
        self.tool = tool
        self.settings = {}


class Plan:
    def __init__(self, steps):
        self.order = [f"n{k}" for k in range(len(steps))]
        self.steps = {f"n{k}": Step(t) for k, t in enumerate(steps)}


def project_cases(projects):
    out = []
    pool = IDS + ["", "unknown", "numpy", "ffmpeg"]
    plans = [[], [None], [Tool(True, [])], [Tool(False, ["rosetta"])],
             [Tool(True, ["ambertools", "pmemd", "ol3", "pmemd", ""]), Tool(True, ["numpy", "cpptraj"])],
             [Tool(True, ["farfar2", "rosetta"]), None, Tool(False, ["dssr"]), Tool(True, ["dssr", "vmd"])]]
    for _ in range(20):
        steps = []
        for _ in range(RNG.randint(1, 6)):
            r = RNG.random()
            if r < 0.1:
                steps.append(None)
            else:
                steps.append(Tool(r < 0.8, [RNG.choice(pool) for _ in range(RNG.randint(0, 6))]))
        plans.append(steps)
    for k, steps in enumerate(plans):
        plain = projects[k % len(projects)]
        project = config.Project.from_dict(plain)
        resources = C.project_resources(project, Plan([t for t in steps]))
        out.append({
            "project": plain,
            "steps": [None if t is None else {"generates": t.generates, "resources": t._resources}
                      for t in steps],
            "generating": [t._resources for t in steps if t is not None and t.generates],
            "resources": resources,
        })
    # no plan: an empty workflow is not resolved
    for plain in projects[:3]:
        project = config.Project.from_dict(plain)
        assert not project.workflow.nodes
        out.append({"project": plain, "steps": [], "generating": [],
                    "resources": C.project_resources(project)})
    return out


def main() -> int:
    out_dir = sys.argv[2]
    os.makedirs(out_dir, exist_ok=True)
    projects = project_variants()
    data = {
        "references": [{"key": r.key, "text": r.text, "bibtex": r.bibtex, "doi": r.doi, "url": r.url,
                        "verified": r.verified} for r in C.REFERENCES.values()],
        "resources": [{"id": r.id, "name": r.name, "keys": list(r.keys), "note": r.note}
                      for r in C.RESOURCES.values()],
        "unique": [{"items": ids, "out": C.unique(ids)} for ids in id_lists()],
        "keys_for": [{"ids": ids, "out": C.keys_for(ids)} for ids in id_lists()],
        "references_of": [{"keys": keys, "out": [r.key for r in C.references(keys)]}
                          for keys in key_lists()],
        "bibtex": [{"keys": keys, "out": C.bibtex(keys)} for keys in key_lists()],
        "text_lines": [{"keys": keys, "out": C.text_lines(keys)} for keys in key_lists()],
        "grouped": [{"ids": ids, "out": [[res.id, [r.key for r in refs]] for res, refs in C.grouped(ids)]}
                    for ids in id_lists()],
        "comment_lines": comment_cases(),
        "pipeline": [{"project": p, "resources": C.pipeline_resources(config.Project.from_dict(p))}
                     for p in projects],
        "project": project_cases(projects),
    }
    path = os.path.join(out_dir, "citation_cases.json")
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        json.dump(data, fh, ensure_ascii=False, indent=1)
        fh.write("\n")
    print(f"wrote {path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
