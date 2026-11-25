"""MDBridge's improvements of project files, external data, the PDB export subset and chains, applied
to MDATK v0.0.4.1 on the test side (docs/DIFFERENCES.md; owner decision D9 of milestone 5; fixes
that later MDATK releases made, named with their version).

The generators run MDATK in-process; where MDBridge deliberately does better, the expected data is
MDATK's output with the improvement applied here, never edited by hand. Each function mirrors the C++
rule exactly (messages, order). Call the apply_* functions after the MDATK checkout is on sys.path;
each is idempotent.

* apply_project()   PRJ-09 / D9: unknown keys of a project file are kept in their table, at any depth,
                    with their values, and to_dict() writes them back after the table's fields
                    (core::projectFromPlain / toPlain); one warning for all of them after a successful
                    from_dict (core::unknownKeysWarning) instead of one per key; --set warns about the
                    unknown keys it adds only (app/commands/overrides.cpp). Patches config._from_plain,
                    config._to_plain, Project.from_dict and, when the command line is importable,
                    cli.apply_override.
* apply_external()  EXT-01, MDATK 0.1.0 (io/external.py:111-121): two external series of one key keep
                    both; the later one gets "-2", "-3" ... (src/io/external.cpp).
* apply_export()    EXP-01, MDATK 0.0.10 (io/pdb_export.py:19-40): an atom subset of commas only is an
                    error; a comma list of plain names with one that is no group reports that group
                    (core::resolveExportSelection).
* apply_chains()    TOP-01, MDATK 0.1.0 (analysis/structure.py:75-84): chains from bonded fragments
                    only when the bonds hold every residue together (core::chainIds).
* apply()           all of them.
"""

from __future__ import annotations

import contextlib
import dataclasses
import re
import typing

_applied = set()


def _once(name):
    if name in _applied:
        return False
    _applied.add(name)
    return True


# ---------------------------------------------------------------------------------------------
# PRJ-09 (D9): unknown keys kept

UNKNOWN = "_mdbridge_unknown_keys"  # the unknown keys of a dataclass object, as read: {key: value}
_BARE = re.compile(r"[A-Za-z0-9_-]+\Z")
_quiet = [0]  # > 0: from_dict logs no warning (cli.apply_override reports new keys itself)
_last = [None]  # the last project from_dict made (apply_override copies its root's unknown keys)


def _key(name, tomlio):
    """core::tomlKey: a bare key when possible, else a basic string."""
    return name if _BARE.match(name) else tomlio.quote(name)


def unknown_key_paths(obj, tomlio, where="project"):
    """core::unknownKeyPaths: a table's own unknown keys, then those inside each field (declaration
    order), depth first; list items as "where[i]"."""
    out = [f"{where}.{_key(k, tomlio)}" for k in getattr(obj, UNKNOWN, {})]
    for f in dataclasses.fields(obj):
        value = getattr(obj, f.name)
        if dataclasses.is_dataclass(value):
            out += unknown_key_paths(value, tomlio, f"{where}.{f.name}")
        elif isinstance(value, list):
            for i, item in enumerate(value):
                if dataclasses.is_dataclass(item):
                    out += unknown_key_paths(item, tomlio, f"{where}.{f.name}[{i}]")
    return out


HINTS = {"project.out_dir": "output.out_dir", "project.cache_dir": "output.cache_dir",
         "project.workflow_dir": "workflow.out_dir"}


def unknown_keys_warning(paths):
    """core::unknownKeysWarning"""
    items = ", ".join(p + (f" (did you mean {HINTS[p]}?)" if p in HINTS else "") for p in paths)
    noun = "unknown key" if len(paths) == 1 else "unknown keys"
    return f"project: {len(paths)} {noun} kept and saved back unchanged (from a newer MDATK, or misspelt?): {items}"


def apply_project():
    if not _once("project"):
        return
    from mdatk import config as C
    from mdatk import tomlio

    def _from_plain(tp, value, where):
        """config._from_plain of v0.0.4.1, unknown keys kept on the object (UNKNOWN)."""
        origin = typing.get_origin(tp)
        args = typing.get_args(tp)
        if origin is typing.Union:
            if value is None:
                return None
            non_none = [a for a in args if a is not type(None)]
            return _from_plain(non_none[0], value, where)
        if dataclasses.is_dataclass(tp):
            if not isinstance(value, dict):
                raise ValueError(f"{where}: expected an object, got {type(value).__name__}")
            hints = typing.get_type_hints(tp)
            names = {f.name for f in dataclasses.fields(tp)}
            kwargs, unknown = {}, {}
            for given, val in value.items():
                key = C._ALIASES.get((tp.__name__, given), given)
                if key not in names:
                    unknown[given] = val
                    continue
                if tp is C.InputParams and key == "trajectories" and isinstance(val, str):
                    val = [val]
                kwargs[key] = _from_plain(hints[key], val, f"{where}.{key}")
            obj = tp(**kwargs)
            if unknown:
                setattr(obj, UNKNOWN, unknown)
            return obj
        return original_from_plain(tp, value, where)

    def _to_plain(obj):
        """config._to_plain of v0.0.4.1, each table's unknown keys after its fields, values as read."""
        if dataclasses.is_dataclass(obj) and not isinstance(obj, type):
            out = {f.name: _to_plain(getattr(obj, f.name)) for f in dataclasses.fields(obj)}
            out.update(getattr(obj, UNKNOWN, {}))
            return out
        if isinstance(obj, (list, tuple)):
            return [_to_plain(v) for v in obj]
        if isinstance(obj, dict):
            return {str(k): _to_plain(v) for k, v in obj.items() if v is not None}
        return obj

    @classmethod
    def from_dict(cls, data, base_dir=None):
        data = dict(data)
        for meta_key in ("format", "format_version", "mdatk_version", "saved", "log", "log_file"):
            data.pop(meta_key, None)
        project = C._from_plain(cls, data, "project")
        paths = unknown_key_paths(project, tomlio)
        if paths and not _quiet[0]:
            C.log.warning("%s", unknown_keys_warning(paths))
        if base_dir:
            project.resolve_paths(base_dir)
        _last[0] = project
        return project

    original_from_plain = C._from_plain
    C._from_plain = _from_plain
    C._to_plain = _to_plain
    C.Project.to_dict = lambda self: C._to_plain(self)
    C.Project.from_dict = from_dict
    try:
        from mdatk import cli
    except ImportError:  # the model generator imports the parameters only
        return
    original_override = cli.apply_override

    def apply_override(project, assignment):
        before = unknown_key_paths(project, tomlio)
        _last[0] = None
        _quiet[0] += 1
        try:
            original_override(project, assignment)
        finally:
            _quiet[0] -= 1
        fresh = _last[0]
        if fresh is not None:  # MDATK copies the fields of the new project; its unknown keys too
            if hasattr(fresh, UNKNOWN):
                setattr(project, UNKNOWN, getattr(fresh, UNKNOWN))
            elif hasattr(project, UNKNOWN):
                delattr(project, UNKNOWN)
        added = [p for p in unknown_key_paths(project, tomlio) if p not in before]
        if added:
            C.log.warning("%s", unknown_keys_warning(added))

    cli.apply_override = apply_override


@contextlib.contextmanager
def quiet_unknown_keys():
    """from_dict without the unknown-keys warning."""
    _quiet[0] += 1
    try:
        yield
    finally:
        _quiet[0] -= 1


# ---------------------------------------------------------------------------------------------
# EXT-01: external series of one key


def apply_external():
    if not _once("external"):
        return
    import os

    import numpy as np
    from mdatk.io import external as E

    def load_external(specs, dataset):
        """io/external.py load_external of MDATK 0.1.0 (its lines 111-121) on v0.0.4.1's readers."""
        result = E.AnalysisResult(name="external", key="", params={}, arrays={}, meta={})
        for n, spec in enumerate(specs):
            if not spec.path:
                continue
            header, data = E.read_table(spec.path, spec.delimiter)
            xi = E._column_index(header, spec.x_column, 0)
            x = data[:, xi]
            unit = spec.x_unit or "ps"
            if unit == "frame":
                target = dataset.frames.astype(float)
            elif unit == "frame1":
                target = dataset.frames.astype(float) + 1.0
            else:
                x = x * E.time_unit_factor(unit)
                target = dataset.time_ps
            order = np.argsort(x, kind="stable")
            x = x[order]
            cols = spec.columns or [h for k, h in enumerate(header) if k != xi]
            label = spec.label or os.path.splitext(os.path.basename(spec.path))[0]
            for col in cols:
                ci = E._column_index(header, col, -1)
                y = data[order, ci]
                values = np.interp(target, x, y, left=np.nan, right=np.nan)
                out = base = E._sanitize(f"{label}.{header[ci]}")
                text = f"{label}: {header[ci]}"
                arr = f"s{n}_{ci}"
                series, k = result.outputs.setdefault("series", {}), 1
                while out in series and series[out]["array"] != arr:
                    k += 1
                    out, text = f"{base}-{k}", f"{label}: {header[ci]} ({k})"
                result.arrays[arr] = values
                result.add_series(out, arr, text, "")
        return result

    E.load_external = load_external
    try:
        from mdatk import pipeline
    except ImportError:
        return
    pipeline.load_external = load_external


# ---------------------------------------------------------------------------------------------
# EXP-01: the PDB export's atom subset


def apply_export():
    if not _once("export"):
        return
    from mdatk.io import pdb_export as X
    from mdatk.selections import SelectionError, split_csv

    def resolve_export_selection(universe, text, groups):
        """io/pdb_export.py of MDATK 0.0.10 (lines 19-40); the unknown group is reported by select()
        as "group:<name>" would be (core::resolveExportSelection)."""
        text = (text or "").strip()
        if not text:
            return universe.atoms
        names = split_csv(text)
        if not names:
            raise SelectionError(f"{text!r} names no residue group (leave the atom subset empty to export "
                                 "every atom)")
        group_names = [g.name for g in groups]
        if all(n in group_names for n in names):
            text = " or ".join(f"group:{n}" for n in names)
        elif "," in text and all(re.fullmatch(r"[A-Za-z0-9_.\-]+", n) for n in names):
            missing = next(n for n in names if n not in group_names)
            return X.select(universe, f"group:{missing}", groups)
        return X.select(universe, text, groups)

    X.resolve_export_selection = resolve_export_selection
    try:
        from mdatk import pipeline
    except ImportError:
        return
    pipeline.resolve_export_selection = resolve_export_selection


# ---------------------------------------------------------------------------------------------
# TOP-01: chains from bonds


def apply_chains():
    if not _once("chains"):
        return
    import numpy as np
    from mdatk.analysis import structure

    def chain_ids(residues):
        """analysis/structure.py chain_ids of MDATK 0.1.0 (lines 54-84)."""
        n = len(residues)
        if n == 0:
            return np.zeros(0, dtype=np.int64)
        first_atoms = residues.atoms[np.unique(residues.atoms.resindices, return_index=True)[1]]
        try:
            cid = np.asarray(first_atoms.chainIDs)
            if len(set(cid)) > 1:
                return np.unique(cid, return_inverse=True)[1].astype(np.int64)
        except (AttributeError, Exception):
            pass
        segs = np.asarray(residues.segids)
        if len(set(segs)) > 1:
            return np.unique(segs, return_inverse=True)[1].astype(np.int64)
        try:
            frag = np.asarray(residues.atoms.fragindices)
        except Exception:
            return np.zeros(n, dtype=np.int64)
        _uniq, first, rows = np.unique(residues.atoms.resindices, return_index=True, return_inverse=True)
        if np.any(frag != frag[first][rows]):
            return np.zeros(n, dtype=np.int64)
        return np.unique(frag[first], return_inverse=True)[1].astype(np.int64)

    structure.chain_ids = chain_ids


def apply():
    apply_project()
    apply_external()
    apply_export()
    apply_chains()
