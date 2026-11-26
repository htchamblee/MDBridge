"""`mdbridge citations`: the specification of MDBridge's new command (milestone 5, owner decision
D14; the format is in docs/specs/m5-plan.md section 7), written on MDATK's own functions so cli_parity.py can record its
expected output from MDATK v0.0.4.1 (with MDBridge's citation improvements, improvements_citations.py,
already applied by the child). MDATK has no such command: its window shows the same text in About and
writes the same BibTeX with *Export citations (BibTeX)…* (mdatk/gui/app.py citation_text,
export_citations); MDBridge adds the run-time citations of the cached results (AP-16, CIT-06).

    mdbridge citations [-p PROJECT] [common options of analyze] [--bibtex FILE] [--json]

* The project: -p (or the defaults) with the overrides of every analysis command (load_project:
  --top, --traj, -o, --cache-dir, --frames, -s ...).
* Keys: the settings' citation keys (keys_for(project_resources(project))) followed by the run-time
  keys of the cached results of the last run (the "citation_keys" of the analyses the cache manifest
  lists, in its order; no manifest or a missing entry: none), without repeats and keys the table
  lacks (core::citations::withRunKeys).
* Without --bibtex / --json: the About text on stdout (MDATK's citation_text, every line followed by
  a newline): "Resources this project uses (N references; each generated script lists only its
  own):", an empty line, per resource with references its name (" - note" when it has one), its
  references indented by four spaces, an empty line; then, when the cached results add references,
  "Recorded by the cached results of the last run:" and those references indented by four spaces
  and an empty line. Without any key: "This project uses no external resource with a reference
  yet." Exit 0.
* --bibtex FILE: writes bibtex(keys) to FILE (UTF-8) and prints "Wrote N references to FILE"; with
  no key nothing is written and "Nothing to cite yet: The project uses no external resource with a
  reference." is printed (exit 0). A file that cannot be written: "error: <OS error>" on stderr,
  exit 1.
* --json: instead of the text, {"keys": [...], "settings_keys": [...], "run_keys": [...] (those of
  the cache, before filtering), "resources": [{"id", "name", "note", "keys"}, ...]} on stdout as JSON
  with an indent of 2 and non-ASCII characters as they are; with --bibtex the file is written first
  and the JSON replaces the "Wrote ..." line.
* Errors of the project (load_project) as every command: "error: <message>" on stderr, exit 1.
"""

from __future__ import annotations

import argparse
import json
import os
import sys

RUN_HEADING = "Recorded by the cached results of the last run:"


def parser(cli):
    p = argparse.ArgumentParser(prog="mdbridge citations")
    cli._common(p)
    p.add_argument("--bibtex", metavar="FILE", help="write the references as BibTeX to FILE")
    p.add_argument("--json", action="store_true", help="print the keys and resources as JSON")
    return p


def cache_run_keys(project):
    """The citation_keys of the analyses of the cache manifest's last run, in its order."""
    from mdatk.io.cache import CacheStore

    try:
        cache = CacheStore(project.cache_dir)
        manifest = cache.read_manifest()
    except (OSError, ValueError):
        return []
    out = []
    for entry in (manifest or {}).get("entries", []):
        loaded = cache.load_result(entry.get("name", ""), entry.get("key", ""))
        if loaded is None:
            return []
        for key in loaded[0].meta.get("citation_keys", []) or []:
            if key not in out:
                out.append(key)
    return out


def text(project_resources, keys, run_only, C):
    if not keys:
        return "This project uses no external resource with a reference yet.\n"
    lines = [f"Resources this project uses ({len(keys)} references; each generated script lists only its own):", ""]
    for res, refs in C.grouped(project_resources):
        if not refs:
            continue
        lines.append(res.name + (f" - {res.note}" if res.note else ""))
        lines += [f"    {r.text}" for r in refs]
        lines.append("")
    if run_only:
        lines.append(RUN_HEADING)
        lines += [f"    {r.text}" for r in C.references(run_only)]
        lines.append("")
    return "\n".join(lines) + "\n"


def main(argv, cli):
    """Runs the command in-process (MDATK imported, its stdout / stderr redirected by the caller)."""
    from mdatk import citations as C
    import improvements_citations

    args = parser(cli).parse_args(argv)
    try:
        project = cli.load_project(args)
        resources = C.project_resources(project)
        settings = C.keys_for(resources)
        run = cache_run_keys(project)
        keys = improvements_citations.with_run_keys(settings, run, C.REFERENCES)
        run_only = [k for k in keys if k not in settings]
        if args.bibtex:
            if not keys:
                print("Nothing to cite yet: The project uses no external resource with a reference.")
                return 0
            with open(args.bibtex, "w", encoding="utf-8") as fh:
                fh.write(C.bibtex(keys))
        if args.json:
            print(json.dumps({"keys": keys, "settings_keys": settings, "run_keys": run,
                              "resources": [{"id": r.id, "name": r.name, "note": r.note, "keys": list(r.keys)}
                                            for r, _refs in C.grouped(resources)]},
                             indent=2, ensure_ascii=False))
        elif args.bibtex:
            print(f"Wrote {len(keys)} references to {args.bibtex}")
        else:
            sys.stdout.write(text(resources, keys, run_only, C))
        return 0
    except (OSError, ValueError, RuntimeError, KeyError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":  # pragma: no cover - for trying it out: MDATK on sys.path
    from mdatk import cli as _cli

    sys.exit(main(sys.argv[1:], _cli))
