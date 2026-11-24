"""What MDATK's data model (mdatk/data.py) computes (tests/data/dataset/dataset_cases.json).

    python3 -I -B make_dataset_cases.py MDATK_CHECKOUT OUT_DIR

Needs numpy; nothing is written into the checkout. Each case builds AnalysisResults (arrays with a
dtype and shape; add_series / add_matrix / add_pairmap calls, with extra keywords; meta entries),
a Dataset of them (frames, time_ps, time_unit) and runs queries on it:

  results[k].outputs_repr     repr(result.meta) after the add_* calls
  queries                     [op, argument, expected]: time / time_label / keys / has / kind_of /
                              series / matrix / pairmap / describe; expected is {"value": ...} or
                              {"error": message} (a KeyError's argument, a ValueError's str())

Also CSRBuilder runs ("csr") and add_* calls that fail ("add_errors"). Floats are float.hex()
strings ("nan", "inf", "-inf"); arrays are {"dtype", "shape", "values"} with values flattened.
"""

import json
import math
import os
import sys

if len(sys.argv) != 3:
    sys.exit("usage: make_dataset_cases.py MDATK_CHECKOUT OUT_DIR")
sys.path.insert(0, os.path.abspath(sys.argv[1]))

import numpy as np  # noqa: E402

from mdatk.data import AnalysisResult, CSRBuilder, Dataset  # noqa: E402


def fhex(x):
    x = float(x)
    if math.isnan(x):
        return "nan"
    if math.isinf(x):
        return "inf" if x > 0 else "-inf"
    return x.hex()


def array_json(a):
    a = np.asarray(a)
    name = a.dtype.name
    flat = a.reshape(-1)
    if name == "bool":
        values = [bool(v) for v in flat]
    elif name.startswith("int"):
        values = [int(v) for v in flat]
    else:
        values = [fhex(v) for v in flat]
    return {"dtype": name, "shape": list(a.shape), "values": values}


def make_array(spec):
    dtype, shape, values = spec
    return np.array(values, dtype=dtype).reshape(shape)


def build(rspec):
    res = AnalysisResult(name=rspec["name"], key=rspec.get("key", ""), params={}, arrays={},
                         meta=dict(rspec.get("meta", {})))
    for name, spec in rspec.get("arrays", []):
        res.arrays[name] = make_array(spec)
    for op in rspec.get("ops", []):
        kind, args, extra = op[0], op[1], op[2] if len(op) > 2 else {}
        getattr(res, kind)(*args, **extra)
    return res


def run(dataset, op, arg):
    if op == "time":
        return [fhex(v) for v in dataset.time(arg)]
    if op == "time_label":
        return dataset.time_label(arg)
    if op == "keys":
        return dataset.keys(arg)
    if op == "has":
        return dataset.has(arg)
    if op == "kind_of":
        return dataset.kind_of(arg)
    if op == "describe":
        return dataset.describe()
    if op == "series":
        s = dataset.series(arg)
        return {"key": s.key, "values": array_json(s.values), "label": s.label, "unit": s.unit,
                "target": None if s.target is None else fhex(s.target), "axis_label": s.axis_label}
    if op == "matrix":
        m = dataset.matrix(arg)
        return {"key": m.key, "values": array_json(m.values), "row_labels": m.row_labels, "label": m.label,
                "kind": m.kind, "classes": m.classes, "unit": m.unit}
    if op == "pairmap":
        p = dataset.pairmap(arg)
        out = {"key": p.key, "offsets": [int(v) for v in p.offsets], "i": [int(v) for v in p.i],
               "j": [int(v) for v in p.j], "values": [fhex(v) for v in p.values], "row_labels": p.row_labels,
               "col_labels": p.col_labels, "label": p.label, "shape": list(p.shape)}
        dense = []
        for f in range(max(p.n_frames, 0)):
            try:
                dense.append({"value": [fhex(v) for v in p.dense(f).reshape(-1)]})
            except IndexError as exc:
                dense.append({"error": str(exc)})
        out["dense"] = dense
        try:
            out["frequency"] = {"value": [fhex(v) for v in p.frequency().reshape(-1)]}
        except IndexError as exc:
            out["frequency"] = {"error": str(exc)}
        return out
    raise ValueError(op)


I32 = "int32"
F32 = "float32"
F64 = "float64"
I64 = "int64"
I8 = "int8"

RMSD = {
    "name": "rmsd", "key": "abc",
    "arrays": [["rmsd", [F64, [4], [0.0, 1.5, float("nan"), 2.25]]],
               ["count", [I32, [4], [1, 2, 3, 4]]],
               ["flags", ["bool", [4], [True, False, True, True]]],
               ["single", [F32, [4], [0.1, 0.2, 1e-8, 3.4e38]]]],
    "ops": [["add_series", ["value", "rmsd", "RMSD", "Å"]],
            ["add_series", ["count", "count", "Contacts", "count"], {"target": 3}],
            ["add_series", ["flags", "flags", "Flag"], {"target": 0.5, "color": "red"}],
            ["add_series", ["single", "single", "Single", ""], {"target": None}],
            ["add_series", ["missing", "nope", "Missing array"]],
            ["add_series", ["value", "rmsd", "RMSD again", "nm"]]],
}
CONTACTS = {
    "name": "contacts",
    "arrays": [["per_res", [I32, [3, 4], list(range(12))]],
               ["occ", [I8, [2, 4], [0, 1, 2, -1, 1, 1, 0, 3]]],
               ["vals", [F64, [2, 4], [0.5] * 8]],
               ["cm_offsets", [I64, [5], [0, 2, 3, 3, 6]]],
               ["cm_i", [I32, [6], [0, 1, 1, 0, 2, 0]]],
               ["cm_j", [I32, [6], [1, 2, 2, 0, -1, 1]]],
               ["cm_v", [F32, [6], [1.0, 0.5, 0.25, 2.0, 1.0, 1.0]]],
               ["bad_offsets", [I64, [3], [0, 1, 2]]],
               ["bad_i", [I32, [2], [5, 0]]],
               ["bad_j", [I32, [2], [0, 0]]],
               ["bad_v", [F32, [2], [1, 1]]],
               ["fl_offsets", [F64, [3], [0.0, 1.0, 2.0]]],
               ["fl_i", [F64, [2], [0.9, 1.0]]],
               ["fl_j", [I64, [2], [0, 1]]],
               ["fl_v", [I32, [2], [3, 4]]],
               ["e_offsets", [I64, [0], []]],
               ["e_i", [I32, [0], []]],
               ["e_j", [I32, [0], []]],
               ["e_v", [F32, [0], []]],
               ["z_offsets", [I64, [1], [0]]],
               ["z_i", [I32, [2], [0, 1]]],
               ["z_j", [I32, [2], [1, 1]]],
               ["z_v", [F32, [2], [1, 1]]]],
    "meta": {"residues": ["A:G1", "A:C2", "B:U3"], "chars": "xyz", "classes_meta": [1, 2]},
    "ops": [["add_matrix", ["per_residue", "per_res", "residues", "Contacts per residue"]],
            ["add_matrix", ["occupancy", "occ", "residues", "Mg occupancy", "class"],
             {"classes": ["none", "inner", "outer", "both"], "unit": ""}],
            ["add_matrix", ["values", "vals", "chars", "Values", "value"], {"unit": "kcal/mol"}],
            ["add_matrix", ["norows", "vals", "absent", "No rows"]],
            ["add_pairmap", ["map", "cm", "residues", "residues", "Contact map"]],
            ["add_pairmap", ["bad", "bad", "residues", "residues", "Bad map"]],
            ["add_pairmap", ["floats", "fl", "residues", "chars", "Float map"], {"note": 1}],
            ["add_pairmap", ["empty", "e", "residues", "residues", "Empty map"]],
            ["add_pairmap", ["zero", "z", "residues", "residues", "Zero frames"]],
            ["add_pairmap", ["noarrays", "none", "residues", "residues", "No arrays"]]],
}
EXTERNAL = {
    "name": "external",
    "arrays": [["s0_1", [F64, [4], [1.0, 2.0, 3.0, 4.0]]]],
    "ops": [["add_series", ["cpp.rmsd", "s0_1", "cpp: rmsd", ""]]],
}
EMPTY = {"name": "empty"}
NO_LABEL = {
    "name": "raw",
    "arrays": [["x", [F64, [4], [4, 3, 2, 1]]]],
    "meta": {"outputs": {"series": {"plain": {"array": "x"}, "none": None}, "matrices": {},
                         "pairmaps": {"m": {"prefix": "p", "rows": "r", "cols": "c"}}}},
}
SLASHED = {
    "name": "a",
    "arrays": [["x", [F64, [4], [1, 1, 1, 1]]]],
    "ops": [["add_series", ["b/c", "x", "Nested"]]],
}

QUERIES = [
    ["time", None], ["time", "ps"], ["time", "ns"], ["time", "us"], ["time", "fs"], ["time", "frame"],
    ["time", "ms"], ["time", ""],
    ["time_label", None], ["time_label", "ps"], ["time_label", "us"], ["time_label", "frame"],
    ["time_label", "fs"], ["time_label", "x"], ["time_label", ""],
    ["keys", "series"], ["keys", "matrices"], ["keys", "pairmaps"], ["keys", "other"],
    ["describe", None],
    ["has", "rmsd/value"], ["has", "rmsd"], ["has", "nope/value"], ["has", "rmsd/nope"],
    ["has", "contacts/map"], ["has", "a/b/c"], ["has", "/"], ["has", ""],
    ["kind_of", "rmsd/value"], ["kind_of", "contacts/per_residue"], ["kind_of", "contacts/map"],
    ["kind_of", "rmsd/nope"], ["kind_of", "x"],
    ["series", "rmsd/value"], ["series", "rmsd/count"], ["series", "rmsd/flags"], ["series", "rmsd/single"],
    ["series", "rmsd/missing"], ["series", "rmsd"], ["series", "nope/x"], ["series", "rmsd/nope"],
    ["series", "contacts/per_residue"], ["series", "external/cpp.rmsd"], ["series", "raw/plain"],
    ["series", "raw/none"], ["series", "a/b/c"], ["series", "empty/x"], ["series", "rmsd/"],
    ["series", "é/x"], ["series", "rmsd/é'\""],
    ["matrix", "contacts/per_residue"], ["matrix", "contacts/occupancy"], ["matrix", "contacts/values"],
    ["matrix", "contacts/norows"], ["matrix", "contacts/map"], ["matrix", "rmsd/value"],
    ["pairmap", "contacts/map"], ["pairmap", "contacts/bad"], ["pairmap", "contacts/floats"],
    ["pairmap", "contacts/empty"], ["pairmap", "contacts/zero"], ["pairmap", "contacts/noarrays"],
    ["pairmap", "raw/m"], ["pairmap", "contacts/per_residue"],
]

CASES = [
    {"name": "full", "frames": [0, 10, 20, 30], "time_ps": [0.0, 100.0, 2000.5, 1e6], "time_unit": "ns",
     "results": [RMSD, CONTACTS, EXTERNAL, EMPTY, NO_LABEL, SLASHED]},
    {"name": "replaced", "frames": [5, 6, 7, 8], "time_ps": [1.0, 2.0, 3.0, 4.0], "time_unit": "frame",
     "results": [RMSD, EXTERNAL, dict(EXTERNAL, arrays=[["s0_1", [F64, [4], [9, 9, 9, 9]]]]), CONTACTS]},
    {"name": "units", "frames": [0, 1], "time_ps": [0.5, 1.5], "time_unit": "us", "results": [EMPTY]},
    {"name": "odd_unit", "frames": [0, 1], "time_ps": [0.5, 1.5], "time_unit": "ms", "results": []},
    {"name": "empty_unit", "frames": [], "time_ps": [], "time_unit": "", "results": []},
]


def plain(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    return value


def cases():
    out = []
    for case in CASES:
        results = [build(r) for r in case["results"]]
        entry = {"name": case["name"], "frames": case["frames"], "time_ps": [fhex(v) for v in case["time_ps"]],
                 "time_unit": case["time_unit"], "results": []}
        for spec, res in zip(case["results"], results):
            entry["results"].append({
                "name": spec["name"], "key": spec.get("key", ""),
                "arrays": [[n, array_json(res.arrays[n])] for n, _ in spec.get("arrays", [])],
                "meta": spec.get("meta", {}),
                "ops": [[op[0], op[1], op[2] if len(op) > 2 else {}] for op in spec.get("ops", [])],
                "outputs_repr": repr(res.meta),
            })
        ds = Dataset(np.array(case["frames"], dtype=np.int64), np.array(case["time_ps"], dtype=np.float64),
                     time_unit=case["time_unit"])
        for res in results:
            ds.add(res)
        queries = []
        for op, arg in QUERIES:
            try:
                queries.append([op, arg, {"value": run(ds, op, arg)}])
            except KeyError as exc:
                queries.append([op, arg, {"error": exc.args[0], "type": "KeyError"}])
            except ValueError as exc:
                queries.append([op, arg, {"error": str(exc), "type": "ValueError"}])
        entry["queries"] = queries
        out.append(entry)
    return out


def csr_cases():
    runs = [
        [],
        [[[], [], []]],
        [[[0, 1], [1, 2], [0.5, 0.25]], [[], [], []], [[2], [0], [1e-8]]],
        [[[0, 2 ** 31], [-1, 3], [0.1, 3.4e38]], [[5, 6, 7], [1, 1, 1], [1, 2, 3]]],
        [[[1, 2, 3], [1], [1.0]]],
    ]
    out = []
    for frames in runs:
        b = CSRBuilder()
        for i, j, v in frames:
            b.add(np.array(i, dtype=np.int64), np.array(j, dtype=np.int64), np.array(v, dtype=np.float64))
        arrays = b.arrays("cm")
        out.append({"frames": [[i, j, [fhex(x) for x in v]] for i, j, v in frames],
                    "arrays": [[k, array_json(a)] for k, a in arrays.items()]})
    return out


def add_error_cases():
    out = []
    for kind, args, extra in [("add_series", ["x", "a", "L"], {"array": "b"}),
                              ("add_series", ["x", "a", "L"], {"unit": "nm"}),
                              ("add_matrix", ["x", "a", "r", "L"], {"rows": "q"}),
                              ("add_pairmap", ["x", "p", "r", "c", "L"], {"prefix": "q"})]:
        res = AnalysisResult("r", "", {}, {})
        try:
            getattr(res, kind)(*args, **extra)
            out.append([kind, args, extra, None])
        except TypeError as exc:
            out.append([kind, args, extra, str(exc)])
    return out


def main():
    out_dir = sys.argv[2]
    os.makedirs(out_dir, exist_ok=True)
    data = {"cases": cases(), "csr": csr_cases(), "add_errors": add_error_cases()}
    with open(os.path.join(out_dir, "dataset_cases.json"), "w", encoding="utf-8") as fh:
        json.dump(data, fh, indent=1, ensure_ascii=False)
        fh.write("\n")


if __name__ == "__main__":
    main()
