#!/usr/bin/env python3
"""The commands of the command line, MDATK v0.0.4.1 against MDBridge.

    python3 -I -B tests/parity/cli_parity.py --mdatk MDATK_CHECKOUT --mdbridge EXE [--only NAME ...]
    python3 -I -B tests/parity/cli_parity.py --mdatk MDATK_CHECKOUT --write tests/data/cli/cli_cases.json
    python3 -I -B tests/parity/cli_parity.py --list

Every case (CASES below) runs `mdatk ARGS` and `EXE ARGS` in two fresh work folders holding the
same inputs (copies of files of tests/data, and files written by the case), then compares the exit
status, stdout, stderr and every file the command left in the folder (byte for byte, except the
analysis cache: see below). `init` / `init-params`, `info`, `list`, `analyze`, `export-pdb` and
`density info` are covered: their results, files (project files, PDB files, the analysis cache),
messages and errors. A case may first run other command lines in the same folder (`before`: an
earlier `analyze` whose cache the case reuses; they must succeed), set environment variables
(`env`: a fake conda installation for `list`) and make files executable (`executable`).

MDATK runs in-process in one child interpreter (this script with --child), from MDATK_CHECKOUT
(nothing is written into it), with MDBridge's independent data tables patched into MDAnalysis
(independent_tables.py, tests/data/README.md): MDBridge deliberately uses its own residue,
element and mass tables, so the reference is MDAnalysis' behaviour with those tables. The child
needs Python 3.13 with MDATK's dependencies at the reference versions (MDAnalysis 2.11.0;
requirements-regenerate.txt); EXE is the `mdbridge` program, the AppImage, or the test build's
`mdbridge_commands_cli` (the same commands without the Qt front end).

Normalised on both sides before comparing (the only differences allowed):
  * the work folder's path -> {work} (each side runs in its own folder);
  * `saved = <time>` of project files -> {time}; run time stamps (2026-10-15_143205) -> {stamp};
  * MDATK's stderr loses the line MDAnalysis logs when netCDF4 is missing (an environment notice,
    not MDATK's output; MDBridge reads NetCDF itself);
  * MDAnalysis' hidden XTC / TRR offset caches (.NAME_offsets.npz / .lock) are not files of the run;
  * MDATK messages that MDBridge words with its own name (NAMING) are adapted;
  * MDBridge's documented improvements (docs/DIFFERENCES.md, improvements_cli.py) are applied to
    MDATK: patched in-process (validation, --frames / --set index messages, the RMSD reference-frame
    warning, external data read before the frame loop, PDB coordinates checked before writing) and
    in its messages (KeyError quotes, stop=None, the out_dir hint, map errors with the path);
  * the analysis cache (ARCHITECTURE: MDBridge keeps its own cache format, mdbridge/analysis/cache.hpp)
    is compared by meaning: file names MDATK_cache.json / MDBridge_cache.json -> {manifest},
    <name>-<key>.npz / <name>-<key>.mdbridge.npz (and .json) -> <name>-{key}.npz (.json); cache
    keys -> {key} in names and messages ("loaded from cache (rmsd-{key})"), since keys hash the
    inputs' mtimes, which differ between the two folders; the manifest and the entries' .json files
    as JSON values with project, version, time stamps, keys and mtime_ns masked ({project},
    "version": {version}, {time}, {key}, 0) and all_entries in a fixed order; the .npz archives as
    their arrays (names in order, dtypes, shapes, values: integers, booleans and strings exactly,
    floats within FLOAT_TOLERANCE relative, at least that much absolute; pair maps <p>_offsets /
    _i / _j / _v per frame as sorted entries - MDBridge lists a frame's stacking pairs sorted, MDATK
    in cKDTree order: tests/analysis/test_mg_stacking.cpp).
Cases marked `exit_only` compare the exit status alone: usage errors (argparse and CLI11 word them
differently; both exit with 2) and documented deviations (DEVIATIONS: map errors that
src/io/include/mdbridge/io/mrc.hpp words differently on purpose). Where MDATK crashes (an exception
its main() does not catch: a Python traceback, exit status 1), MDBridge must exit with 1 and print
"error: <message>"; stdout and files are compared as usual.

Help texts are not compared (argparse and CLI11 lay them out differently; MDBridge's `--help` lists
only the commands it has).

--write runs MDATK alone and writes the normalised results as the expected data of
tests/cli/test_commands.cpp (tests/data/cli/cli_cases.json); tests/parity/regenerate_all.py runs it.
Files up to 16 KiB are stored as text, larger ones as their SHA-256 and size; cache files as
{"json": masked value} and {"npz": [[name, dtype, shape, flat values as one line of JSON text], ...]}
(floats as numbers, non-finite ones as "nan" / "inf" / "-inf"); several entries of one analysis
with different keys as {"several": [...]}, matched in any order. Standard library only (the child
imports MDATK; .npz archives are read with zipfile). Exit status 1 on any difference.
"""

from __future__ import annotations

import argparse
import ast
import difflib
import hashlib
import io
import json
import math
import os
import re
import shutil
import stat
import struct
import subprocess
import sys
import tempfile
import traceback
import zipfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)  # -I leaves the script's folder out
import improvements_cli  # noqa: E402

DATA = os.path.normpath(os.path.join(HERE, "..", "data"))
LIMIT = 16384

DEMO = {"demo_dry.parm7": "topology/demo_dry.parm7", "demo_dry.nc": "io/trajectories/demo_dry.nc"}
DEMO_RUN = ["--top", "demo_dry.parm7", "--traj", "demo_dry.nc"]
NO_STAMP = ["-s", "output.timestamp_files=false"]
MAPS = "io/maps/files"

# Maps whose error MDBridge words differently on purpose beyond what improvements_cli.py adjusts:
# exit status only. (Broken gzip members raise MapError "{path}: <gzip's message>" where Python's
# gzip.BadGzipFile has no path: compared after adding the path, CLI-04. Maps on which MDATK crashes
# are compared as described above.)
DEVIATIONS = set()

PROJECT = """format = "mdatk-project"
notes = "relative paths, resolved against this file's folder"

[inputs]
topology = "demo_dry.parm7"
trajectories = ["demo_dry.nc"]
dt_ps = 5.0
time_unit = "ps"

[[groups]]
name = "Aptamer"
residues = "1-4, 21-24"
color = "#123456"

[[groups]]
name = "Target"
selection = "resid 5:20 and not name H*"

[output]
out_dir = "results"
timestamp_files = false

[output.pdb]
frames = "1, 3"
times = "100"
selection = "Aptamer"
prefix = "snap"
"""

GROUPS = json.dumps([
    {"name": "A", "residues": "1-3"},
    {"name": "bad name", "selection": "resname ZZZ"},
    {"name": "A", "residues": "20-22"},
    {"name": "Ions", "selection": "resname MG"},
    {"name": "Broken", "selection": "name ("},
    {"name": "Chain", "residues": "X:5-6"},
    {"name": "Empty", "residues": ""},
    {"name": "Ünï", "residues": "7"},
])


def case(name, argv, inputs=None, setup=None, exit_only=False, before=None, env=None, executable=None):
    """inputs: work-folder path -> file of tests/data copied there; setup: path -> text written there
    ({work}: the folder); before: command lines run first in the folder; env: environment variables
    of the runs ({work}: the folder; "" counts as unset for MDATK); executable: setup files made
    executable."""
    return {"name": name, "argv": argv, "inputs": dict(inputs or {}), "setup": dict(setup or {}),
            "exit_only": exit_only, "before": list(before or []), "env": dict(env or {}),
            "executable": list(executable or [])}


# `list` without conda installations, and with a fake one: an installation whose base and three
# environments provide cpptraj, one environment without it and one whose cpptraj is not executable,
# a mamba root with a non-ASCII environment name, and environments registered in
# ~/.conda/environments.txt (one also found under its root, one missing). Every variable MDATK reads
# is set, so the machine's own installations stay out (except any in COMMON_ROOTS' absolute folders).
NO_CONDA = {"HOME": "{work}/home", "CONDA_EXE": "", "CONDA_PREFIX": "", "MAMBA_ROOT_PREFIX": "", "CONDA_ROOT": ""}
FAKE_CONDA = {**NO_CONDA, "CONDA_EXE": "{work}/miniforge/condabin/conda",
              "CONDA_PREFIX": "{work}/miniforge/envs/py312", "MAMBA_ROOT_PREFIX": "{work}/mamba"}
CONDA_EXECUTABLES = ["miniforge/bin/cpptraj", "miniforge/envs/AmberTools26/bin/cpptraj",
                     "miniforge/envs/amber24/bin/cpptraj", "miniforge/envs/zeta/bin/cpptraj",
                     "mamba/envs/Ünï/bin/cpptraj", "other/envs/ambertools/bin/cpptraj",
                     "standalone/bin/cpptraj"]
CONDA_LAYOUT = {
    **{path: "#!/bin/sh\n" for path in CONDA_EXECUTABLES},
    "miniforge/condabin/conda": "#!/bin/sh\n",
    "miniforge/envs/py312/conda-meta/history": "",
    "miniforge/envs/noexec/bin/cpptraj": "not executable\n",
    "home/.conda/environments.txt": "{work}/other/envs/ambertools\n\n  {work}/standalone  \r\n{work}/missing\n"
                                    "{work}/miniforge/envs/zeta\n",
}
ALL_GROUP = ["-s", 'groups=[{"name": "All", "selection": "all"}]', "-s", "analyses.contacts.enabled=false",
             "-s", "analyses.hbonds.enabled=false", "-s", "analyses.mg.enabled=false", "-s",
             "analyses.stacking.enabled=false", "-s", "analyses.energy.enabled=false"]
EXTERNAL = "time_ps,potential,temp K\n0,1.5,300\n100,2.5,301\n200,4.0,299.5\n300,3.0,300.25\n"


def map_cases():
    out = []
    for name in sorted(os.listdir(os.path.join(DATA, MAPS))):
        out.append(case(f"density_info_{name}", ["density", "info", f"maps/{name}"],
                        {f"maps/{name}": f"{MAPS}/{name}"}, exit_only=name in DEVIATIONS))
    return out


CASES = [
    # ---- init / init-params
    case("init_default", ["init", "new.mdatk.toml"]),
    case("init_alias_options", ["init-params", "p.mdatk.toml", *DEMO_RUN, "--frames", "2:20:3", "-o", "out",
                                "--cache-dir", "cache", "--cpptraj", "/opt/amber/bin/cpptraj", "--conda-env",
                                "AmberTools26", "--clean", "-s", "analyses.contacts.cutoff=5", "-s",
                                "groups.0.residues=1-6", "-s", "notes=5", "-s", 'output.charts.formats=["png","svg"]',
                                "-s", "inputs.dt_ps=2.5", "-s", "groups.-1.name=Last", "-s",
                                'layout.panels.3.options={"mode": "mean"}', "-s", "visual.rotation.1=90",
                                "-s", "cleaning.enabled=false", "--no-clean", "--clean"], DEMO),
    case("init_relative", ["init", "rel.mdatk.toml", "--top", "sub/../demo_dry.parm7", "--traj", "a.nc", "b.nc",
                           "--frames", ":7", "-s", "inputs.trajectories=[\"x.nc\"]", "-s", "output.out_dir=o"]),
    case("init_from_project", ["init", "copy.mdatk.toml", "-p", "base.mdatk.toml", "--no-clean"], DEMO,
         {"base.mdatk.toml": PROJECT}),
    case("init_from_params_log", ["init", "replay.mdatk.toml", "--params", "MDATK_params_log.toml"],
         {"MDATK_params_log.toml": "projects/MDATK_params_log.toml"}),
    case("init_exists", ["init", "exists.mdatk.toml"], setup={"exists.mdatk.toml": "keep me\n"}),
    case("init_overwrite", ["init", "base.mdatk.toml", "--overwrite", "-s", "notes=rewritten"],
         setup={"base.mdatk.toml": PROJECT}),
    case("init_json", ["init", "preset.json"]),
    case("init_no_folder", ["init", "missing/new.mdatk.toml"]),
    case("init_set_no_equals", ["init", "x.mdatk.toml", "-s", "nonsense"]),
    case("init_set_empty_key", ["init", "x.mdatk.toml", "-s", " . =3"]),
    case("init_set_unknown", ["init", "x.mdatk.toml", "-s", "analyses.bogus=1"]),
    case("init_set_free_form", ["init", "x.mdatk.toml", "-s", "layout.panels.3.options.mode=x"]),
    case("init_set_bool", ["init", "x.mdatk.toml", "-s", "cleaning.enabled=1"]),
    case("init_set_index", ["init", "x.mdatk.toml", "-s", "groups.7.name=x"]),
    case("init_set_index_text", ["init", "x.mdatk.toml", "-s", "groups.first.name=x"]),
    case("init_set_index_dashes", ["init", "x.mdatk.toml", "-s", "groups.--1.name=x"]),
    case("init_set_type", ["init", "x.mdatk.toml", "-s", "analyses.contacts.cutoff=far"]),
    case("init_set_table", ["init", "x.mdatk.toml", "-s", 'groups={"name": "x"}']),
    case("init_set_property", ["init", "x.mdatk.toml", "-s", "out_dir=5"]),
    case("init_set_method", ["init", "x.mdatk.toml", "-s", "validate.x=5"]),
    case("init_set_section", ["init", "x.mdatk.toml", "-s", 'inputs={"topology": "t.parm7", "bogus": 1}']),
    case("init_frames_bad", ["init", "x.mdatk.toml", "--frames", "1:x"]),
    case("init_frames_bad_part", ["init", "x.mdatk.toml", "--frames", "1::2x", "-s", "groups.--1.name=x"]),
    case("init_set_cache_dir", ["init", "x.mdatk.toml", "-s", "cache_dir=c", "-s", "workflow_dir=w"]),
    # a single-value option given twice: the last value (argparse); --traj: the last list
    case("init_options_twice", ["init", "x.mdatk.toml", "--top", "a.parm7", "--top", "b.parm7", "--traj", "a.nc",
                                "b.nc", "--traj", "c.nc", "--frames", "1:2", "--frames", "3:9:2", "-o", "o1", "-o",
                                "o2", "--cache-dir", "c1", "--cache-dir", "c2", "--cpptraj", "x", "--cpptraj", "y",
                                "--conda-env", "e1", "--conda-env", "e2", "-p", "missing.mdatk.toml", "--params",
                                "base.mdatk.toml"], DEMO, {"base.mdatk.toml": PROJECT}),
    case("init_missing_output", ["init"], exit_only=True),
    case("init_unknown_option", ["init", "x.mdatk.toml", "--bogus"], exit_only=True),

    # ---- info
    case("info_demo", ["info", *DEMO_RUN], DEMO),
    case("info_no_topology", ["info"]),
    case("info_groups", ["info", *DEMO_RUN, "-s", f"groups={GROUPS}"], DEMO),
    case("info_project", ["info", "-p", "base.mdatk.toml"], DEMO, {"base.mdatk.toml": PROJECT}),
    case("info_two_trajectories", ["info", "--top", "demo_dry.parm7", "--traj", "demo_dry.nc", "demo_dry.xtc",
                                   "-s", "inputs.dt_ps=4"],
         {**DEMO, "demo_dry.xtc": "io/trajectories/files/demo_dry.xtc"}),
    case("info_dcd", ["info", "--top", "demo_dry.parm7", "--traj", "demo_dry.dcd", "--frames", "::0"],
         {**DEMO, "demo_dry.dcd": "io/trajectories/files/demo_dry.dcd"}),
    case("info_pdb", ["info", "--top", "elements.pdb", "-s", 'groups=[{"name": "All", "selection": "all"}]'],
         {"elements.pdb": "io/loader/files/elements.pdb"}),
    case("info_gro", ["info", "--top", "precise.gro"], {"precise.gro": "io/loader/files/precise.gro"}),
    case("info_missing_topology", ["info", "--top", "missing.parm7"]),
    case("info_missing_trajectory", ["info", "--top", "demo_dry.parm7", "--traj", "missing.nc"], DEMO),
    case("info_parm7_alone", ["info", "--top", "demo_dry.parm7"], DEMO),
    case("info_invalid_parameters", ["info", *DEMO_RUN, "-s", "analyses.hbonds.distance_cutoff=-1", "-s",
                                     "analyses.hbonds.angle_cutoff=181", "-s", "analyses.mg.outer_cutoff=0", "-s",
                                     "analyses.mg.inner_cutoff=3", "-s", "analyses.energy.dielectric=0", "-s",
                                     "analyses.stacking.planarity_cutoff=-0.5", "-s",
                                     "analyses.proximity_search_radius=0", "-s",
                                     'analyses.rmsd.targets=[{"label": "A B"}, {"label": "C"}, {"label": "A_B"}, '
                                     '{"label": "!"}, {"label": "target3"}]'], DEMO),
    case("info_unknown_option", ["info", "--bogus"], exit_only=True),

    # ---- export-pdb
    case("export_demo", ["export-pdb", *DEMO_RUN, "--frame-list", "0, 10, 20-23:2, last", "--selection", "Aptamer",
                         *NO_STAMP], DEMO),
    case("export_times", ["export-pdb", *DEMO_RUN, "--times", "55ps, 0.1ns, 10ns", "--selection", "resname MG",
                          "-o", "pdbs", *NO_STAMP, "-s", "output.pdb.prefix=ion"], DEMO),
    case("export_stamp", ["export-pdb", *DEMO_RUN, "--frame-list", "3", "--selection", "Target"], DEMO),
    case("export_all_atoms", ["export-pdb", *DEMO_RUN, "--frame-list", "last", *NO_STAMP], DEMO),
    case("export_project", ["export-pdb", "-p", "base.mdatk.toml"], DEMO, {"base.mdatk.toml": PROJECT}),
    case("export_project_override", ["export-pdb", "-p", "base.mdatk.toml", "--frame-list", "", "--times",
                                     "0.01ns", "--selection", "Aptamer, Target", "-q"], DEMO,
         {"base.mdatk.toml": PROJECT}),
    case("export_verbose_layout", ["export-pdb", *DEMO_RUN, "--frame-list", "2", "-v", *NO_STAMP, "-s",
                                   "layout.panels.1.type=structure", "--selection", "resid 3"], DEMO),
    case("export_no_frames", ["export-pdb", *DEMO_RUN, *NO_STAMP], DEMO),
    case("export_invalid", ["export-pdb", "--top", "demo_dry.parm7", "--frame-list", "0"], DEMO),
    case("export_bad_selection", ["export-pdb", *DEMO_RUN, "--frame-list", "0", "--selection", "name ("], DEMO),
    case("export_empty_selection", ["export-pdb", *DEMO_RUN, "--frame-list", "0", "--selection", "resname ZZZ"],
         DEMO),
    case("export_bad_frames", ["export-pdb", *DEMO_RUN, "--frame-list", "0, 99"], DEMO),
    case("export_options_twice", ["export-pdb", *DEMO_RUN, "--frame-list", "0", "--frame-list", "2, 4", "--selection",
                                  "Target", "--selection", "Aptamer", *NO_STAMP], DEMO),
    # a frame beyond the PDB coordinate limits: nothing is written (CLI-08)
    case("export_coordinate_limits", ["export-pdb", "--top", "limit_models.pdb", "--traj", "limit_models.pdb",
                                      "--frame-list", "0, 1, 2, 3", *NO_STAMP, *ALL_GROUP],
         {"limit_models.pdb": "io/pdb/export/files/limit_models.pdb"}),
    case("export_pdb_trajectory", ["export-pdb", "--top", "chains.pdb", "--traj", "chains.pdb", "--frame-list",
                                   "0", *NO_STAMP, "-s", 'groups=[{"name": "All", "selection": "all"}]',
                                   "-s", "analyses.contacts.enabled=false", "-s", "analyses.hbonds.enabled=false",
                                   "-s", "analyses.mg.enabled=false", "-s", "analyses.stacking.enabled=false",
                                   "-s", "analyses.energy.enabled=false"],
         {"chains.pdb": "io/pdb/export/files/chains.pdb"}),

    # ---- list
    case("list", ["list"], env=NO_CONDA),
    case("list_conda", ["list"], setup=CONDA_LAYOUT, env=FAKE_CONDA, executable=CONDA_EXECUTABLES),
    case("list_extra_argument", ["list", "x"], exit_only=True),

    # ---- analyze
    case("analyze_demo", ["analyze", *DEMO_RUN], DEMO),
    case("analyze_only", ["analyze", *DEMO_RUN, "--only",
                          " rmsd, contacts@Target+Aptamer,rmsd ,stacking@Target+Aptamer "], DEMO),
    case("analyze_only_internal", ["analyze", *DEMO_RUN, "--only", "trace,proximity,atoms,basepairs", "-q",
                                   "--frames", "5:6"], DEMO),
    case("analyze_only_twice", ["analyze", *DEMO_RUN, "--only", "mg", "--force", "--only", "rmsd", "--force"], DEMO),
    case("analyze_only_empty", ["analyze", *DEMO_RUN, "--only", "", "-s", "output.video.enabled=false"], DEMO),
    case("analyze_reuse", ["analyze", *DEMO_RUN, "--frames", "::3"], DEMO,
         before=[["analyze", *DEMO_RUN, "--frames", "::3", "-q"]]),
    case("analyze_reuse_partial", ["analyze", *DEMO_RUN, "--only", "rmsd,hbonds,energy"], DEMO,
         before=[["analyze", *DEMO_RUN, "--only", "rmsd,mg", "-q"]]),
    case("analyze_reuse_changed", ["analyze", *DEMO_RUN, "--only", "contacts,rmsd", "-s",
                                   "analyses.rmsd.enabled=false"], DEMO,
         before=[["analyze", *DEMO_RUN, "--only", "rmsd,contacts", "-s", "analyses.contacts.cutoff=4", "-q"]]),
    case("analyze_force", ["analyze", *DEMO_RUN, "--force", "--only", "rmsd,stacking"], DEMO,
         before=[["analyze", *DEMO_RUN, "--only", "rmsd,stacking", "-q"]]),
    case("analyze_frames", ["analyze", *DEMO_RUN, "--frames", "2:20:3", "--only", "rmsd,contacts,trace"], DEMO),
    case("analyze_frames_negative", ["analyze", *DEMO_RUN, "--frames=-5::2", "-s", "inputs.dt_ps=2.5",
                                     "--only", "hbonds"], DEMO),
    case("analyze_frames_empty", ["analyze", *DEMO_RUN, "--frames", "30:"], DEMO),
    case("analyze_dirs", ["analyze", *DEMO_RUN, "-o", "out", "--cache-dir", "cache", "--only", "mg"], DEMO),
    case("analyze_use_pbc", ["analyze", *DEMO_RUN, "-s", "analyses.use_pbc=true", "--only", "contacts,mg,stacking"],
         DEMO),
    case("analyze_layout", ["analyze", *DEMO_RUN, "-v", "--frames", "3:4", "-s", "layout.panels.1.type=structure",
                            "-s", 'visual.representations=[{"style": "licorice", "groups": ["Aptamer"]}]', "-s",
                            'output.charts.series=["rmsd@Aptamer/Backbone", "external/x", "bogus@Nope/x"]',
                            "-s", "analyses.contacts.enabled=false", "-s", "analyses.energy.enabled=false"], DEMO),
    case("analyze_nothing", ["analyze", *DEMO_RUN, "-s", "output.video.enabled=false", "-s",
                             "output.charts.enabled=false",
                             *(arg for a in ("contacts", "energy", "hbonds", "mg", "rmsd", "stacking")
                               for arg in ("-s", f"analyses.{a}.enabled=false"))], DEMO),
    case("analyze_external", ["analyze", *DEMO_RUN, "--only", "rmsd", "-s",
                              'inputs.external=[{"path": "ext.csv", "label": "MD"}]'], DEMO, {"ext.csv": EXTERNAL}),
    case("analyze_external_bad", ["analyze", *DEMO_RUN, "--only", "rmsd", "-s",
                                  'inputs.external=[{"path": "bad.csv"}]'], DEMO, {"bad.csv": "a,b\n1,x\n"}),
    case("analyze_external_missing", ["analyze", *DEMO_RUN, "--only", "rmsd", "-s",
                                      'inputs.external=[{"path": "missing.csv"}]'], DEMO),
    case("analyze_project", ["analyze", "-p", "base.mdatk.toml", "--only", "rmsd@Aptamer,hbonds"], DEMO,
         {"base.mdatk.toml": PROJECT}),
    case("analyze_unknown", ["analyze", *DEMO_RUN, "--only", "rmsd,nope"], DEMO),
    case("analyze_unknown_group", ["analyze", *DEMO_RUN, "--only", "contacts@Nope"], DEMO),
    case("analyze_invalid", ["analyze", *DEMO_RUN, "-s", "analyses.contacts.cutoff=-1", "-s",
                             "layout.preset=bogus"], DEMO),
    case("analyze_no_topology", ["analyze"]),
    case("analyze_missing_trajectory", ["analyze", "--top", "demo_dry.parm7", "--traj", "missing.nc"], DEMO),
    case("analyze_no_trajectory", ["analyze", "--top", "demo_dry.parm7", "--only", "rmsd"], DEMO),
    case("analyze_options_twice", ["analyze", "--top", "missing.parm7", "--top", "demo_dry.parm7", "--traj",
                                   "missing.nc", "--traj", "demo_dry.nc", "--frames", "0:3", "--frames", "2:8:2",
                                   "--only", "rmsd", "-o", "o1", "-o", "out", "--cache-dir", "c1", "--cache-dir",
                                   "cache"], DEMO),
    case("analyze_reference_beyond", ["analyze", *DEMO_RUN, "--only", "rmsd", "--frames", "::4", "-s",
                                      "analyses.rmsd.reference_frame=999"], DEMO),
    case("analyze_reference_negative", ["analyze", *DEMO_RUN, "--only", "rmsd", "--frames", "::6", "-s",
                                        "analyses.rmsd.reference_frame=-3"], DEMO),
    case("analyze_invalid_parameters", ["analyze", *DEMO_RUN, "-s", "analyses.stacking.angle_cutoff=-1", "-s",
                                        "analyses.energy.cutoff=0", "-s", "analyses.mg.angle_tolerance=0", "-s",
                                        'analyses.rmsd.targets=[{"label": "x"}, {"label": "x"}]'], DEMO),
    case("analyze_unknown_option", ["analyze", "--bogus"], exit_only=True),

    # ---- density info
    *map_cases(),
    case("density_missing", ["density", "info", "missing.mrc"]),
    case("density_folder", ["density", "info", "."]),
    case("density_no_action", ["density"], exit_only=True),
]

STAMP_RE = re.compile(r"\d{4}-\d{2}-\d{2}_\d{6}")
SAVED_RE = re.compile(r"^saved = \d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}(\.\d+)?$", re.M)
OFFSETS_RE = re.compile(r"^\..+_offsets\.(npz|lock)$")
NETCDF_NOTICE = "WARNING: netCDF4 is not available. Writing AMBER ncdf files will be slow.\n"
TRACEBACK = "Traceback (most recent call last):"
# MDATK's words that MDBridge's libraries replace by its own name (applied to MDATK's output)
NAMING = [("JSON presets are no longer supported; MDATK 0.0.4.1 uses TOML project files",
           "JSON presets are no longer supported; MDBridge uses TOML project files"),
          ("Run 'mdatk analyze' to list the keys of a system.", "Run 'mdbridge analyze' to list the keys of a system.")]
# the analysis cache: keys in names and messages, file names, the float tolerance of .npz arrays
CACHE_KEY_RE = re.compile(r"-[0-9a-f]{16}(?=\)|(\.mdbridge)?\.(npz|json)\b)")
MANIFEST_RE = re.compile(r"\b(MDATK|MDBridge)_cache\.json\b")
CACHE_ENTRY_RE = re.compile(r"-\{key\}\.mdbridge\.(npz|json)$")
FLOAT_TOLERANCE = 1e-9
NPY_TYPES = {"<f8": "d", "<f4": "f", "<i8": "q", "<i4": "i", "|i1": "b", "|b1": "?"}


# ---------------------------------------------------------------------------------------------
# work folders

def with_work(text, folder):
    return text.replace("{work}", os.path.abspath(folder))


def prepare(case_, folder):
    os.makedirs(folder)
    for dest, source in case_["inputs"].items():
        target = os.path.join(folder, dest)
        os.makedirs(os.path.dirname(target), exist_ok=True)
        shutil.copyfile(os.path.join(DATA, source), target)
    for dest, text in case_["setup"].items():
        target = os.path.join(folder, dest)
        os.makedirs(os.path.dirname(target), exist_ok=True)
        with open(target, "w", encoding="utf-8", newline="\n") as fh:
            fh.write(with_work(text, folder))
    for dest in case_["executable"]:
        target = os.path.join(folder, dest)
        os.chmod(target, os.stat(target).st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)


def case_env(case_, folder):
    return {k: with_work(v, folder) for k, v in case_["env"].items()}


def normalise(text, folder):
    for path in {os.path.realpath(folder), os.path.abspath(folder)}:
        text = text.replace(path, "{work}")
    return MANIFEST_RE.sub("{manifest}", CACHE_KEY_RE.sub("-{key}", STAMP_RE.sub("{stamp}", text)))


def normalise_name(rel, folder):
    """A file's relative path as compared: {work}, {stamp}, the cache's names."""
    rel = normalise(rel, folder)
    return CACHE_ENTRY_RE.sub(r"-{key}.\1", rel)


def collect(case_, folder):
    """Files left in the folder (inputs included when the run changed them): compared path -> [bytes]
    (several for the cache entries of one analysis with different keys)."""
    out = {}
    for here, _dirs, names in os.walk(folder):
        for name in names:
            if OFFSETS_RE.match(name):
                continue
            full = os.path.join(here, name)
            rel = os.path.relpath(full, folder).replace(os.sep, "/")
            with open(full, "rb") as fh:
                data = fh.read()
            if rel in case_["inputs"]:
                continue  # copies of test data: never written by these commands
            if rel in case_["setup"] and data == with_work(case_["setup"][rel], folder).encode("utf-8"):
                continue
            key = normalise_name(rel, folder)
            if key in out and not is_cache_file(key):
                raise RuntimeError(f"{case_['name']}: {rel} and another file are both compared as {key}")
            out.setdefault(key, []).append(data)  # cache entries of one analysis with several keys
    return dict(sorted(out.items()))


def is_cache_file(rel):
    return rel.endswith("{manifest}") or rel.endswith("-{key}.npz") or rel.endswith("-{key}.json")


def masked_cache_json(value, folder, top=True):
    """A cache manifest or entry description with what legitimately differs masked."""
    if isinstance(value, dict):
        out = {}
        for k, v in value.items():
            if top and k in ("mdatk_version", "mdbridge_version"):
                out["version"] = "{version}"
            elif top and k == "project":
                out[k] = "{project}"
            elif top and k in ("updated", "created"):
                out[k] = "{time}"
            elif k == "key" and isinstance(v, str):
                out[k] = "{key}"
            elif k == "mtime_ns":
                out[k] = 0
            elif k == "file" and isinstance(v, str):
                out[k] = normalise_name(v, folder)
            else:
                out[k] = masked_cache_json(v, folder, False)
        if "all_entries" in out and isinstance(out["all_entries"], list):
            out["all_entries"] = sorted(out["all_entries"], key=lambda e: json.dumps(e, sort_keys=True))
        return out
    if isinstance(value, list):
        return [masked_cache_json(v, folder, False) for v in value]
    if isinstance(value, str):
        return normalise(value, folder)
    return value


def json_float(x):
    return x if math.isfinite(x) else ("nan" if math.isnan(x) else ("inf" if x > 0 else "-inf"))


def read_npz(data):
    """[[name, dtype, shape, flat values], ...] of a .npz archive (C order arrays of NPY_TYPES / <U)."""
    out = []
    with zipfile.ZipFile(io.BytesIO(data)) as archive:
        for member in archive.namelist():
            raw = archive.read(member)
            assert raw[:6] == b"\x93NUMPY", member
            if raw[6] == 1:
                (length,), start = struct.unpack("<H", raw[8:10]), 10
            else:
                (length,), start = struct.unpack("<I", raw[8:12]), 12
            header = ast.literal_eval(raw[start:start + length].decode("latin1"))
            assert not header["fortran_order"], member
            body = raw[start + length:]
            dtype, shape = header["descr"], list(header["shape"])
            count = math.prod(shape)
            if dtype.startswith("<U"):
                width = int(dtype[2:])
                values = [body[i * 4 * width:(i + 1) * 4 * width].decode("utf-32-le").rstrip("\0")
                          for i in range(count)]
            else:
                values = list(struct.unpack(f"<{count}{NPY_TYPES[dtype]}", body[:count * struct.calcsize(
                    NPY_TYPES[dtype])]))
                if dtype in ("<f8", "<f4"):
                    values = [json_float(v) for v in values]
            out.append([member[:-4] if member.endswith(".npy") else member, dtype, shape, values])
    return out


def file_entry(rel, data, folder):
    if is_cache_file(rel):
        if rel.endswith(".npz"):
            # the flat values as one line of JSON text (compact; the C++ test parses it)
            return {"npz": [[name, dtype, shape, json.dumps(values, ensure_ascii=False, separators=(",", ":"))]
                            for name, dtype, shape, values in read_npz(data)]}
        return {"json": masked_cache_json(json.loads(data.decode("utf-8")), folder)}
    try:
        text = data.decode("utf-8")
    except UnicodeDecodeError:
        return {"sha256": hashlib.sha256(data).hexdigest(), "size": len(data)}
    text = SAVED_RE.sub("saved = {time}", normalise(text, folder))
    if len(text.encode("utf-8")) <= LIMIT:
        return {"text": text}
    raw = text.encode("utf-8")
    return {"sha256": hashlib.sha256(raw).hexdigest(), "size": len(raw)}


def written_entry(rel, datas, folder):
    """file_entry, or {"several": [file_entry, ...]} (in a fixed order) for several cache entries."""
    if len(datas) == 1:
        return file_entry(rel, datas[0], folder)
    return {"several": sorted((file_entry(rel, d, folder) for d in datas), key=lambda e: json.dumps(e, sort_keys=True))}


def result(case_, folder, code, stdout, stderr, before):
    files = collect(case_, folder)
    # results that hold absolute paths or OS error texts, or read a project file's paths: defined on
    # Linux only (the C++ test skips them on Windows). The cache manifest's path that analyze prints
    # and the paths in cache files are compared on Windows too (the C++ test normalises their form).
    texts = [stderr, *(data[0].decode("utf-8", "replace") for rel, data in files.items() if not is_cache_file(rel))]
    if case_["argv"][:1] != ["analyze"]:
        texts.append(stdout)
    return {"exit": code, "stdout": normalise(stdout, folder), "stderr": normalise(stderr, folder),
            "written": {rel: written_entry(rel, data, folder) for rel, data in files.items()},
            "before_exits": before,
            "posix_only": any("{work}" in normalise(t, folder) for t in texts) or "[Errno" in stderr
            or any(a in ("-p", "--project", "--params") for a in case_["argv"]) or bool(case_["env"])}


# ---------------------------------------------------------------------------------------------
# MDATK, in-process in a child interpreter

def child(mdatk, jobs_file, results_file):
    import logging
    import warnings

    warnings.filterwarnings("ignore")
    sys.path.insert(0, os.path.abspath(mdatk))
    sys.path.insert(0, HERE)
    import independent_tables

    independent_tables.apply()
    import improvements_analyses

    improvements_analyses.apply_mdatk()  # the analysis improvements (docs/DIFFERENCES.md)

    import improvements_cli

    improvements_cli.apply()  # MDBridge's documented improvements (docs/DIFFERENCES.md)
    from mdatk import cli

    def run(argv):
        out, err = io.StringIO(), io.StringIO()
        root = logging.getLogger()
        for handler in list(root.handlers):  # basicConfig() in cli.main binds the current sys.stderr
            root.removeHandler(handler)
        old = sys.stdout, sys.stderr
        sys.stdout, sys.stderr = out, err
        try:
            code = cli.main(argv)
        except SystemExit as exc:  # argparse
            code = exc.code if isinstance(exc.code, int) else 1
        except BaseException:  # noqa: BLE001 - an uncaught exception: Python prints it, exit 1
            traceback.print_exc()
            code = 1
        finally:
            sys.stdout, sys.stderr = old
        return [code, out.getvalue(), err.getvalue().replace(NETCDF_NOTICE, "")]

    with open(jobs_file, encoding="utf-8") as fh:
        jobs = json.load(fh)
    results = []
    start = os.getcwd()
    for argv, folder, before, env in jobs:
        saved = dict(os.environ)
        os.environ.update(env)
        os.chdir(folder)
        try:
            exits = [run(b)[0] for b in before]
            results.append([*run(argv), exits])
        finally:
            os.chdir(start)
            os.environ.clear()
            os.environ.update(saved)
    with open(results_file, "w", encoding="utf-8") as fh:
        json.dump(results, fh)


def run_mdatk(cases, root, mdatk, python):
    jobs = []
    for c in cases:
        folder = os.path.join(root, "mdatk", c["name"])
        prepare(c, folder)
        jobs.append([c["argv"], folder, c["before"], case_env(c, folder)])
    jobs_file, results_file = os.path.join(root, "jobs.json"), os.path.join(root, "results.json")
    with open(jobs_file, "w", encoding="utf-8") as fh:
        json.dump(jobs, fh)
    subprocess.run([python, "-I", "-B", os.path.abspath(__file__), "--child", mdatk, jobs_file, results_file],
                   check=True)
    with open(results_file, encoding="utf-8") as fh:
        outputs = json.load(fh)
    results = {}
    for c, job, (code, stdout, stderr, before) in zip(cases, jobs, outputs):
        if any(before):
            raise RuntimeError(f"{c['name']}: a command line run before the case failed in MDATK ({before})")
        for old, new in NAMING:
            stdout, stderr = stdout.replace(old, new), stderr.replace(old, new)
        stdout, stderr = improvements_cli.improve_text(stdout), improvements_cli.improve_text(stderr)
        stderr = improved_map_error(c, stderr)
        res = result(c, job[1], code, stdout, stderr, before)
        if TRACEBACK in stderr:
            # MDATK crashed (an exception its main() does not catch: MemoryError, zlib.error,
            # AttributeError ...): Python prints the traceback and exits with 1, MDBridge reports
            # "error: <message>"; the exit status, stdout and files are compared
            res["stderr"] = None
            res["mdatk_crash"] = stderr.strip().splitlines()[-1]
        results[c["name"]] = res
    return results


def improved_map_error(case_, stderr):
    """`density info`: MapError texts as mrc.cpp words them (CLI-04: always the path, byte-order
    problems named as such)."""
    if case_["argv"][:2] != ["density", "info"] or not stderr.startswith("error: ") or stderr.count("\n") != 1:
        return stderr
    path = case_["argv"][2]
    source = case_["inputs"].get(path)
    if source is None:
        return stderr
    with open(os.path.join(DATA, source), "rb") as fh:
        message = improvements_cli.improve_map_message(stderr[len("error: "):-1], fh.read(), path)
    if not message.startswith(path + ": ") and not message.startswith("[Errno"):
        message = f"{path}: {message}"
    return f"error: {message}\n"


def run_mdbridge(cases, root, exe):
    out = {}
    for c in cases:
        folder = os.path.join(root, "mdbridge", c["name"])
        prepare(c, folder)
        env = dict(os.environ, APPIMAGE_EXTRACT_AND_RUN="1", **case_env(c, folder))
        before = [subprocess.run([*exe, *b], cwd=folder, capture_output=True, env=env, timeout=300).returncode
                  for b in c["before"]]
        run = subprocess.run([*exe, *c["argv"]], cwd=folder, capture_output=True, env=env, timeout=300)
        decode = lambda b: b.decode("utf-8", "surrogateescape")  # noqa: E731
        out[c["name"]] = result(c, folder, run.returncode, decode(run.stdout), decode(run.stderr), before)
    return out


# ---------------------------------------------------------------------------------------------

def close(a, b):
    if isinstance(a, str) or isinstance(b, str):  # non-finite values
        return a == b
    return abs(a - b) <= FLOAT_TOLERANCE * max(1.0, abs(a), abs(b))


def npz_differences(rel, expected, actual):
    """Arrays compared as described above; [] when they agree."""
    names = [a[0] for a in expected]
    if names != [a[0] for a in actual]:
        return [f"{rel}: arrays {names} / {[a[0] for a in actual]}"]
    e_arrays = {n: (n, t, s, json.loads(v)) for n, t, s, v in expected}
    a_arrays = {n: (n, t, s, json.loads(v)) for n, t, s, v in actual}
    problems = []
    pairmaps = {n[:-len("_offsets")] for n in names if n.endswith("_offsets")
                and all(n[:-len("_offsets")] + s in e_arrays for s in ("_i", "_j", "_v"))}
    for name in names:
        (_n, e_type, e_shape, e_values), (_m, a_type, a_shape, a_values) = e_arrays[name], a_arrays[name]
        if (e_type, e_shape) != (a_type, a_shape):
            problems.append(f"{rel}: {name} is {e_type} {e_shape} / {a_type} {a_shape}")
        elif any(name == p + s for p in pairmaps for s in ("_i", "_j", "_v")):
            continue  # compared with the offsets
        elif name[:-len("_offsets")] in pairmaps and name.endswith("_offsets"):
            p = name[:-len("_offsets")]
            if e_values != a_values:
                problems.append(f"{rel}: {name} differs")
                continue
            for f in range(len(e_values) - 1):
                lo, hi = e_values[f], e_values[f + 1]
                rows = [sorted(zip(*(arrays[p + s][3][lo:hi] for s in ("_i", "_j", "_v"))))
                        for arrays in (e_arrays, a_arrays)]
                if any(x[:2] != y[:2] or not close(x[2], y[2]) for x, y in zip(*rows)):
                    problems.append(f"{rel}: pair map {p}, frame {f}: {rows[0][:6]}... / {rows[1][:6]}...")
                    break
        elif e_type in ("<f8", "<f4"):
            bad = [i for i, (x, y) in enumerate(zip(e_values, a_values)) if not close(x, y)]
            if bad:
                problems.append(f"{rel}: {name}[{bad[0]}] {e_values[bad[0]]} / {a_values[bad[0]]} "
                                f"({len(bad)} values beyond the tolerance)")
        elif e_values != a_values:
            bad = next(i for i, (x, y) in enumerate(zip(e_values, a_values)) if x != y)
            problems.append(f"{rel}: {name}[{bad}] {e_values[bad]!r} / {a_values[bad]!r}")
    return problems


def differences(c, expected, actual):
    if expected["exit"] != actual["exit"]:
        problems = [f"exit status: MDATK {expected['exit']}, MDBridge {actual['exit']}"]
        if c["exit_only"]:
            return problems
    else:
        problems = []
    if c["exit_only"]:
        return []
    if expected["before_exits"] != actual["before_exits"]:
        problems.append(f"command lines run first: exit statuses {expected['before_exits']} / "
                        f"{actual['before_exits']}")
    if expected["stderr"] is None and not actual["stderr"].startswith("error: "):
        problems.append(f"MDATK crashed ({expected['mdatk_crash']}), MDBridge printed no error: "
                        f"{actual['stderr']!r}")
    for part in ("stdout", "stderr"):
        if expected[part] is not None and expected[part] != actual[part]:
            problems.append(f"{part} differs\n{diff(expected[part], actual[part])}")
    for rel in sorted(set(expected["written"]) | set(actual["written"])):
        e, a = expected["written"].get(rel), actual["written"].get(rel)
        if e is None or a is None:
            problems.append(f"{rel}: written by {'MDBridge' if e is None else 'MDATK'} only")
            continue
        if "several" in e or "several" in a:
            es, unused = e.get("several", [e]), list(a.get("several", [a]))
            for entry in es:  # each matched with one of MDBridge's, in any order
                match = next((x for x in unused if not entry_differences(rel, entry, x)), None)
                if match is None:
                    problems.append(f"{rel}: MDATK's entry has no counterpart in MDBridge's {len(unused)}")
                else:
                    unused.remove(match)
            if len(es) != len(a.get("several", [a])):
                problems.append(f"{rel}: {len(es)} files in MDATK, {len(a.get('several', [a]))} in MDBridge")
        else:
            problems.extend(entry_differences(rel, e, a))
    return problems


def entry_differences(rel, e, a):
    if "npz" in e and "npz" in a:
        return npz_differences(rel, e["npz"], a["npz"])
    if "json" in e and "json" in a:
        if e["json"] == a["json"]:
            return []
        return [f"{rel} differs\n" + diff(json.dumps(e["json"], indent=1, ensure_ascii=False),
                                          json.dumps(a["json"], indent=1, ensure_ascii=False))]
    if e == a:
        return []
    if "text" in e and "text" in a:
        return [f"{rel} differs\n{diff(e['text'], a['text'])}"]
    return [f"{rel} differs: {e} / {a}"]


def diff(expected, actual):
    lines = difflib.unified_diff(expected.splitlines(), actual.splitlines(), "MDATK", "MDBridge", lineterm="", n=2)
    return "\n".join(list(lines)[:40])


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--mdatk", help="checkout of MDATK v0.0.4.1")
    parser.add_argument("--python", default=sys.executable, help="interpreter that runs MDATK (default: this one)")
    parser.add_argument("--mdbridge", help="the mdbridge program, AppImage or mdbridge_commands_cli to compare")
    parser.add_argument("--write", help="write MDATK's normalised results to this JSON file instead")
    parser.add_argument("--only", nargs="+", help="run these cases only")
    parser.add_argument("--keep", help="run in this (new) folder and keep it")
    parser.add_argument("--list", action="store_true", help="list the cases")
    parser.add_argument("--child", nargs=3, help=argparse.SUPPRESS)
    args = parser.parse_args(argv)
    if args.child:
        child(*args.child)
        return 0
    cases = [c for c in CASES if not args.only or c["name"] in args.only]
    if args.list:
        for c in cases:
            print(f"{c['name']:40s} {' '.join(c['argv'])}")
        return 0
    if not args.mdatk or bool(args.mdbridge) == bool(args.write):
        parser.error("give --mdatk and one of --mdbridge / --write")
    root = os.path.abspath(args.keep) if args.keep else tempfile.mkdtemp(prefix="cli_parity_")
    os.makedirs(root, exist_ok=True)
    try:
        expected = run_mdatk(cases, root, args.mdatk, args.python)
        if args.write:
            data = {"cases": [{**c, **expected[c["name"]]} for c in cases]}
            os.makedirs(os.path.dirname(os.path.abspath(args.write)), exist_ok=True)
            with open(args.write, "w", encoding="utf-8", newline="\n") as fh:
                json.dump(data, fh, indent=1, ensure_ascii=False)
                fh.write("\n")
            print(f"{len(cases)} cases written to {args.write}")
            return 0
        # the program runs in each case's folder: a path relative to here would be looked up there
        exe = os.path.abspath(args.mdbridge) if os.path.dirname(args.mdbridge) else args.mdbridge
        actual = run_mdbridge(cases, root, [exe])
        failed = 0
        for c in cases:
            problems = differences(c, expected[c["name"]], actual[c["name"]])
            print(f"{'ok' if not problems else 'FAIL':4s} {c['name']}")
            for problem in problems:
                print("     " + problem.replace("\n", "\n     "))
            failed += bool(problems)
        print(f"{len(cases) - failed} of {len(cases)} cases give the same results in MDATK and MDBridge")
        return 1 if failed else 0
    finally:
        if not args.keep:
            shutil.rmtree(root, ignore_errors=True)


if __name__ == "__main__":
    sys.exit(main())
