"""MDBridge's improvements of project validation, parameters and the command line, applied to MDATK
v0.0.4.1 on the test side (docs/IMPROVEMENTS.md sections 3 and 5, docs/DIFFERENCES.md).

The generators and cli_parity.py run MDATK in-process; where MDBridge deliberately does better, the
expected data is MDATK's output with the improvement applied here, never edited by hand. Pure Python
(numpy only for the PDB coordinate check). Two kinds of adjustment:

* apply(): patches MDATK's functions in-process (behaviour that changes what else happens):
    PRJ-01 / PRJ-03  Project.validate: the analysis parameter checks of MDBridge's
                     core::analysisParameterProblems take the place of MDATK's contacts / Mg checks
    RN-20            Project.validate: visual.renderer "vulkan" is accepted (owner decision D7) and
                     the unknown-renderer message lists it: "(auto, software, vtk or vulkan)."
    PRJ-02           RMSD setup warns when reference_frame is clamped
    PRJ-04           apply_override: "groups.--1.name: no list element '--1'"
    PRJ-05           load_project: "--frames: 'x' is not a frame number or range"
    SEL-17           info selects the groups with periodic distance keywords only with use_pbc
                     (MDATK 0.1.2)
    CLI-02           Pipeline.analyze reads the external data once before the frame loop
    CLI-08           export_frames checks every frame's coordinates before writing any file
    CMD-01           build_demo checks the frame count, the time step and the sequence first
    CMD-08           build_script: the cpptraj script's first line names MDBridge and its version
* improve_text() / improve_strings(): messages reworded in place:
    CLI-01           a KeyError printed without the quotes of str(KeyError)
    CLI-03           "frame selection start=30 stop=None ..." -> "stop=end"
  and improve_map_message() (CLI-04) for the MRC reader's errors.
(PRJ-07's hint is part of the unknown-keys warning of improvements_core_io.py, PRJ-09.)

Each function mirrors the C++ rule exactly (messages, order); see the item ids in the C++ sources.
"""

from __future__ import annotations

import ast
import copy
import gzip
import io
import math
import os
import re
import struct

# ---------------------------------------------------------------------------------------------
# PRJ-01 / PRJ-03: validate()


def rmsd_data_key(label, k):
    """_out_key of MDATK's RMSD analysis."""
    return re.sub(r"[^A-Za-z0-9_.+\-]+", "_", label).strip("_") or f"target{k}"


def analysis_parameter_problems(project):
    """core::analysisParameterProblems (src/core/project.cpp), in its order and words."""
    out = []
    a = project.analyses

    def positive(value, what):
        if not value > 0:
            out.append(f"{what} must be positive.")
            return False
        if not math.isfinite(value):
            out.append(f"{what} must be finite.")
            return False
        return True

    def angle(value, what):
        if not 0 <= value <= 180:
            out.append(f"{what} must be between 0 and 180 degrees.")

    positive(a.contacts.cutoff, "Contact cutoff")
    positive(a.hbonds.distance_cutoff, "H-bond distance cutoff")
    angle(a.hbonds.angle_cutoff, "H-bond angle cutoff")
    inner = positive(a.mg.inner_cutoff, "Mg inner-sphere cutoff")
    outer = positive(a.mg.outer_cutoff, "Mg outer-sphere cutoff")
    if inner and outer and a.mg.inner_cutoff > a.mg.outer_cutoff:
        out.append("Mg inner-sphere cutoff must not exceed the outer-sphere cutoff.")
    positive(a.mg.angle_tolerance, "Mg angle tolerance")
    positive(a.stacking.distance_cutoff, "Stacking distance cutoff")
    angle(a.stacking.angle_cutoff, "Stacking angle cutoff")
    positive(a.stacking.offset_cutoff, "Stacking offset cutoff")
    positive(a.stacking.planarity_cutoff, "Stacking planarity cutoff")
    positive(a.energy.cutoff, "Energy cutoff")
    positive(a.energy.dielectric, "Energy dielectric constant")
    positive(a.proximity_search_radius, "Proximity search radius")
    targets = a.rmsd.targets
    for j in range(1, len(targets)):
        key = rmsd_data_key(targets[j].label, j)
        for i in range(j):
            if rmsd_data_key(targets[i].label, i) == key:
                out.append(f"RMSD targets {i + 1} and {j + 1} give the same data key {key!r}; rename one.")
                break
    return out


def improved_validate(original):
    """MDATK's validate() with its two analysis checks (contacts cutoff <= 0, Mg inner > outer)
    replaced, at their place, by analysis_parameter_problems()."""

    def validate(self):
        clean = copy.deepcopy(self)  # neither of MDATK's two checks fires
        clean.analyses.contacts.cutoff, clean.analyses.mg.inner_cutoff, clean.analyses.mg.outer_cutoff = 1.0, 1.0, 2.0
        problems = original(clean)
        probe = copy.deepcopy(clean)
        probe.analyses.contacts.cutoff = -1.0
        at = original(probe).index("Contact cutoff must be positive.")
        return problems[:at] + analysis_parameter_problems(self) + problems[at:]

    return validate


def renderer_validate(original):
    """RN-20: MDATK's validate() with "vulkan" a valid visual.renderer and the message naming it."""
    old_message = "(auto, software or vtk)."

    def validate(self):
        if self.visual.renderer == "vulkan":
            clean = copy.deepcopy(self)
            clean.visual.renderer = "auto"
            return original(clean)
        return [p.replace(old_message, "(auto, software, vtk or vulkan).") if p.startswith("Unknown renderer ")
                else p for p in original(self)]

    return validate


# ---------------------------------------------------------------------------------------------
# messages (CLI-01, CLI-03)

FRAMES_RE = re.compile(r"(frame selection start=-?\d+ stop=)None( step=)")
QUOTED_RE = re.compile(r"^(error: |WARNING: skipping data source '[^'\n]*': )(['\"].*['\"])$", re.M)


def _unquote(match):
    try:
        value = ast.literal_eval(match.group(2))
    except (ValueError, SyntaxError):
        return match.group(0)
    return match.group(1) + value if isinstance(value, str) else match.group(0)


def improve_text(text):
    """MDATK's console / log text with the reworded messages of CLI-01 and CLI-03."""
    if not isinstance(text, str):
        return text
    text = FRAMES_RE.sub(r"\1end\2", text)
    return QUOTED_RE.sub(_unquote, text)


def improve_strings(value):
    """improve_text on every string of a JSON-like value."""
    if isinstance(value, str):
        return improve_text(value)
    if isinstance(value, list):
        return [improve_strings(v) for v in value]
    if isinstance(value, dict):
        return {k: improve_strings(v) for k, v in value.items()}
    return value


# ---------------------------------------------------------------------------------------------
# CLI-04: map errors

MODES = {0, 1, 2, 6, 12}
BYTE_ORDER = "cannot determine the byte order of the map header"
MODE_RE = re.compile(r"^(.*): unsupported MRC mode (-?\d+) \(supported: 0, 1, 2, 6, 12\)$")
ANGLES_RE = re.compile(r"^(.*): non-orthogonal cell angles (\[[^\]]*\]) are not supported$")


def _header(data):
    if data[:2] == b"\x1f\x8b":
        try:
            return gzip.GzipFile(fileobj=io.BytesIO(data)).read(1024)
        except (OSError, EOFError):
            return b""
    return data[:1024]


def _byte_order(raw):
    """'<' or '>' as mrc.cpp chooses it (stamp, else the order whose MODE and MAPC look valid)."""
    if raw[212] == 0x44:
        return "<"
    if raw[212] == 0x11:
        return ">"
    for order in ("<", ">"):
        mode = struct.unpack(order + "i", raw[12:16])[0]
        crs = sorted(struct.unpack(order + "3i", raw[64:76]))
        if mode in MODES and crs in ([1, 2, 3], [0, 0, 0]):
            return order
    return None


def improve_map_message(message, data, path):
    """The MapError text of mrc.cpp for MDATK's `message` about the map `data` (bytes) named `path`
    in messages: the byte-order message carries the path, and a MODE or cell angles that make sense
    only in the other byte order are reported as a byte-order problem."""
    if message == BYTE_ORDER:
        return f"{path}: {message}"
    raw = _header(data)
    if len(raw) < 1024:
        return message
    order = _byte_order(raw)
    if order is None:
        return message
    other = ">" if order == "<" else "<"
    m = MODE_RE.match(message)
    if m and struct.unpack(other + "i", raw[12:16])[0] in MODES:
        return f"{m.group(1)}: {BYTE_ORDER} (mode {m.group(2)})"
    m = ANGLES_RE.match(message)
    if m and all(abs(a - 90.0) <= 0.01 for a in struct.unpack(other + "3f", raw[52:64])):
        return f"{m.group(1)}: {BYTE_ORDER} (cell angles {m.group(2)})"
    return message


# ---------------------------------------------------------------------------------------------
# CLI-08: PDB export


def check_export_coordinates(universe, frames, atoms):
    """pdb_export.cpp's check before writing (exportFrames): the first requested frame with a
    coordinate outside the PDB limits raises ValueError."""
    import numpy as np

    low, high = np.float32(-999.9995), np.float32(9999.9995)
    for f in frames:
        universe.trajectory[int(f)]
        pos = atoms.positions.astype(np.float32)
        if not np.all((low < pos) & (pos <= high)):
            raise ValueError("PDB files must have coordinate values between -1000.000 and 9999.999 "
                             f"Angstroem: frame {int(f)} is outside this range; no files were written.")


# ---------------------------------------------------------------------------------------------
# apply()

_applied = False


def apply():
    """Patch the imported MDATK modules (idempotent). Call after the MDATK checkout is on sys.path."""
    global _applied
    if _applied:
        return
    _applied = True
    import improvements_core_io  # project files, external data, export subset, chains (fix wave B)

    improvements_core_io.apply()
    from mdatk import config

    config.Project.validate = renderer_validate(improved_validate(config.Project.validate))

    try:
        from mdatk import cli
    except ImportError:  # the model generator imports the parameters only
        return
    apply_rmsd()
    _patch_cli(cli)
    _patch_pipeline()
    _patch_export()
    _patch_demo()
    _patch_cpptraj_script()


def apply_rmsd():
    """Only the RMSD reference-frame warning (PRJ-02), for the analysis generators."""
    global _rmsd_applied
    if not _rmsd_applied:
        _rmsd_applied = True
        _patch_rmsd()


_rmsd_applied = False


def _patch_rmsd():  # PRJ-02
    from mdatk.analysis import rmsd

    original = rmsd.RMSDAnalysis.setup

    def setup(self, frame_indices):
        p = self.project.analyses.rmsd
        n = len(frame_indices)
        if p.targets and n > 0:
            if p.reference_frame >= n:
                self.ctx.warn(f"RMSD: reference frame {p.reference_frame} is beyond the {n} analysed frames; "
                              "using the last")
            elif p.reference_frame < 0:
                self.ctx.warn(f"RMSD: reference frame {p.reference_frame} is negative; using the first")
        return original(self, frame_indices)

    rmsd.RMSDAnalysis.setup = setup


def _patch_cli(cli):
    original_info = cli.cmd_info

    def cmd_info(args):  # SEL-17 (MDATK 0.1.2): the groups are selected periodic only with use_pbc
        import improvements_selections
        from mdatk import selections

        load, select = cli.load_project, selections.select
        loaded = {}

        def load_project(a):
            loaded["project"] = load(a)
            return loaded["project"]

        def select_open(universe, text, groups=(), allow_empty=False):
            periodic = bool(loaded["project"].analyses.use_pbc)
            return improvements_selections.select_mdbridge(selections, universe, text, groups, allow_empty,
                                                           periodic=periodic)

        cli.load_project, selections.select = load_project, select_open
        try:
            return original_info(args)
        finally:
            cli.load_project, selections.select = load, select

    cli.cmd_info = cmd_info
    original_override = cli.apply_override

    def apply_override(project, assignment):  # PRJ-04
        try:
            return original_override(project, assignment)
        except ValueError as exc:
            m = re.match(r"invalid literal for int\(\) with base 10: (.*)$", str(exc))
            if not m or "=" not in assignment:
                raise
            path = assignment.split("=", 1)[0]
            part = ast.literal_eval(m.group(1))
            raise ValueError(f"{path}: no list element {part!r}") from None

    original_load = cli.load_project

    def load_project(args, resolve=True):  # PRJ-05
        try:
            return original_load(args, resolve)
        except ValueError:
            frames = getattr(args, "frames", None)
            if frames:
                for part in (frames.split(":") + ["", "", ""])[:3]:
                    try:
                        int(part) if part else 0
                    except ValueError:
                        raise ValueError(f"--frames: {frames!r} is not a frame number or range") from None
            raise

    cli.apply_override = apply_override
    cli.load_project = load_project


def _patch_pipeline():  # CLI-02
    from mdatk import pipeline
    from mdatk.data import Dataset
    import numpy as np

    original = pipeline.run_analyses

    def run_analyses(ctx, objs, frame_indices, reporter=None, *args, **kwargs):
        p = ctx.project
        if p.inputs.external:
            empty = Dataset(np.array([], dtype=np.int64), np.array([], dtype=np.float64), {},
                            time_unit=p.inputs.time_unit)
            try:
                pipeline.load_external(p.inputs.external, empty)
            except (OSError, ValueError) as exc:
                raise pipeline.PipelineError(f"external data: {exc}") from exc
        return original(ctx, objs, frame_indices, reporter, *args, **kwargs)

    pipeline.run_analyses = run_analyses


def _patch_export():  # CLI-08
    from mdatk import pipeline
    from mdatk.io import pdb_export

    original = pdb_export.export_frames

    def export_frames(universe, frames, atoms, out_dir, prefix, progress=None, stamp=""):
        if frames:
            os.makedirs(out_dir, exist_ok=True)
            check_export_coordinates(universe, frames, atoms)
        return original(universe, frames, atoms, out_dir, prefix, progress=progress, stamp=stamp)

    pdb_export.export_frames = export_frames
    pipeline.export_frames = export_frames


def _patch_demo():  # CMD-01 (src/app/commands/demo.cpp buildDemo)
    from mdatk import demo

    original = demo.build_demo

    def build_demo(out_dir, n_frames=120, dt_ps=10.0, solvent=True, sequence="GCAUGGCAUCGA", seed=7):
        if n_frames < 1:
            raise ValueError(f"the demo needs at least one frame (--n-frames {n_frames})")
        if not (math.isfinite(dt_ps) and dt_ps > 0):
            raise ValueError(f"the time between frames must be a positive number of ps (--dt {dt_ps!r})")
        if len(sequence) < 4 or any(c not in "ACGU" for c in sequence):
            raise ValueError(f"demo sequence {sequence!r}: at least four of the letters A, C, G and U")
        return original(out_dir, n_frames, dt_ps, solvent, sequence, seed)

    demo.build_demo = build_demo


def _patch_cpptraj_script():  # CMD-08 (src/analysis/cleaning.cpp buildCpptrajScript)
    from mdatk.cleaning import cpptraj

    original = cpptraj.build_script
    header = f"# Generated by {cpptraj.PROJECT_ABBREV} {cpptraj.__version__} on "

    def build_script(params, topology, trajectories, out_dir):
        script = original(params, topology, trajectories, out_dir)
        assert script.startswith(header), script.splitlines()[0]
        # the version as the normalisation writes it (cli_parity.py GENERATED_RE: MDBridge <version>
        # -> MDBridge {version}); the generators run from a copy of tests/ and tools/ only
        return "# Generated by MDBridge {version} on " + script[len(header):]

    cpptraj.build_script = build_script  # clean() looks it up in its module

