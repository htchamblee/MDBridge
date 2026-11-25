"""What MDATK's commands `render`, `charts` and `preview` print and write (mdatk/cli.py cmd_render,
cmd_charts, cmd_preview; mdatk/pipeline.py), for MDBridge's commands (src/app/commands/render.cpp;
docs/specs/m4-plan.md sections 1.7, 1.8 and 3).

    /opt/mdbridge-venv/bin/python -I -B make_render_cli_cases.py MDATK_CHECKOUT OUT_DIR

Needs MDATK's dependencies and FFmpeg's `ffmpeg` and `ffprobe` on PATH (the regenerate job
installs Ubuntu 24.04's package); nothing is written into the checkout. MDATK runs in-process
(mdatk.cli.main, as cli_parity.py's child runs it), patched by render_common.setup() (independent
tables, analysis / CLI / rendering improvements), each case in a fresh work folder holding copies of
tests/data/topology/demo_dry.parm7 and tests/data/io/trajectories/demo_dry.nc and the project file
demo.mdatk.toml (written by MDATK's Project.save, then its `saved` time fixed). Writes
OUT_DIR/cli/render_cli_cases.json:

  project_file  the text of demo.mdatk.toml every case starts from (the demo project on the demo
                files, out_dir "results", 320 x 180 video at dpi 40 and 12 fps, charts at dpi 40 of
                three series, time-stamped file names; the cases change it with --set)
  cases         per case: name, description, before (command lines run first in the same folder,
                e.g. `analyze` for --from-cache; they must succeed), argv, and
                exit          the exit status (1 for an uncaught exception: MDATK crashed)
                stdout        normalised (below); the commands print JSON or a path
                stderr        normalised: log lines, progress lines and warnings; for a crash only
                              the traceback's last line ("crash")
                files         every file left in the folder but the inputs, normalised path ->
                              info: size_px (PNG), size_pt (SVG width / height, PDF MediaBox),
                              magic (first bytes), and for .mp4 videos what ffprobe reports
                              (stream: codec_type, codec_name, width, height, pix_fmt,
                              r_frame_rate, avg_frame_rate, time_base, nb_read_frames; format:
                              format_name, nb_streams; audio streams) plus faststart (the moov
                              box before mdat); for MDATK_params_log.toml its run's `outputs`
                              table and `timestamp`; the analysis cache only by name
                step_log      the project file's step log entries after the run (MDATK's
                              projectfile.read_log; time fields masked)

Normalisation (as cli_parity.py): the work folder -> {work}; run time stamps (2026-10-15_143205) ->
{stamp}; cache keys -> {key}; MDATK_cache.json -> {manifest}; "Finished in N s" -> "Finished in
{seconds} s"; ISO times in the step log -> {time}, durations -> 0.0. The video's bytes are not
recorded (encoders differ; MDBridge's encoder and codec are its own, RN-07); its codec name is
MDATK's (libx264 from Ubuntu's FFmpeg) and is not compared by the C++ tests.
"""

import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import render_common as rc  # noqa: E402  (before numpy: pins its kernels)

CHECKOUT, OUT = rc.args("make_render_cli_cases.py MDATK_CHECKOUT OUT_DIR")
rc.setup(CHECKOUT, collect=False)

import io  # noqa: E402
import json  # noqa: E402
import logging  # noqa: E402
import re  # noqa: E402
import shutil  # noqa: E402
import struct  # noqa: E402
import subprocess  # noqa: E402
import tempfile  # noqa: E402
import tomllib  # noqa: E402
import traceback  # noqa: E402

import improvements_render  # noqa: E402
from mdatk import cli  # noqa: E402
from mdatk.config import PanelSpec  # noqa: E402
from mdatk.demo import demo_project  # noqa: E402
from mdatk.projectfile import read_log  # noqa: E402

PROJECT = "demo.mdatk.toml"
INPUTS = {"demo_dry.parm7": os.path.join(rc.DATA, "topology", "demo_dry.parm7"),
          "demo_dry.nc": os.path.join(rc.DATA, "io", "trajectories", "demo_dry.nc")}
P = ["-p", PROJECT]
STAMP_RE = re.compile(r"\d{4}-\d{2}-\d{2}_\d{6}")
CACHE_KEY_RE = re.compile(r"-[0-9a-f]{16}(?=\)|(\.mdbridge)?\.(npz|json)\b)")
MANIFEST_RE = re.compile(r"\b(MDATK|MDBridge)_cache\.json\b")
FINISHED_RE = re.compile(r"Finished in [0-9.]+ s")
SAVED_RE = re.compile(r"^saved = .*$", re.M)
ISO_RE = re.compile(r"^\d{4}-\d{2}-\d{2}[T ]\d{2}:\d{2}:\d{2}")
NETCDF_NOTICE = "WARNING: netCDF4 is not available. Writing AMBER ncdf files will be slow.\n"
TRACEBACK = "Traceback (most recent call last):"


def s(key, value):
    return ["-s", f"{key}={json.dumps(value)}"]


CASES = [
    # (name, before, argv, description)
    ("render", [], ["render", *P], "full run: analyses, video, charts (time-stamped names), params log, step log"),
    ("render_no_stamp", [], ["render", *P, "--no-video", *s("output.timestamp_files", False)],
     "charts only, plain names"),
    ("render_from_cache", [["analyze", *P]], ["render", *P, "--from-cache"],
     "video and charts from the cache manifest"),
    ("render_from_cache_no_video", [["analyze", *P]], ["render", *P, "--from-cache", "--no-video"],
     "charts from the cache"),
    ("render_from_cache_nothing", [["analyze", *P]], ["render", *P, "--from-cache", "--no-video", "--no-charts"],
     "nothing enabled: {}"),
    ("render_stride", [], ["render", *P, "--no-charts", *s("output.video.frame_stride", 5)],
     "every 5th frame: 0, 5, ..., 20 (the last frame is skipped, RN-17)"),
    ("render_mpeg4", [], ["render", *P, "--no-charts", *s("output.video.codec", "mpeg4"),
                          *s("output.video.filename", "clip.mp4"), *s("output.video.fps", 30)],
     "codec mpeg4, own file name, 30 fps"),
    ("render_codec_fallback", [], ["render", *P, "--no-charts", *s("output.video.codec", "nosuchcodec")],
     "an encoder FFmpeg lacks: warning and fallback (RN-07)"),
    ("render_odd_size", [], ["render", *P, "--no-charts", *s("output.video.width", 333),
                             *s("output.video.height", 187), *s("output.video.dpi", 37)],
     "odd size and a DPI that does not divide it (RN-06)"),
    ("render_no_ffmpeg", [], ["render", *P, "--no-charts", *s("output.video.ffmpeg_path", "/nonexistent/ffmpeg")],
     "ffmpeg missing: VideoError"),
    ("render_pdb", [], ["render", *P, "--no-video", "--no-charts", "--pdb", *s("output.pdb.frames", "0, last")],
     "PDB export only"),
    ("charts", [], ["charts", *P], "charts after the analyses"),
    ("charts_svg", [["analyze", *P]], ["charts", *P, "--from-cache", "--formats", "svg,pdf", "--dpi", "30"],
     "SVG and PDF from the cache"),
    ("charts_series", [["analyze", *P]], ["charts", *P, "--from-cache", "--series", "energy/total,nokey/x"],
     "own series list with a missing key"),
    ("charts_bad_format", [["analyze", *P]], ["charts", *P, "--from-cache", "--formats", "png,bogus"],
     "an unknown format: Matplotlib's ValueError after the first file (RN-11)"),
    ("preview", [], ["preview", *P], "frame 0 at the video size: results/MDATK_preview_00000.png"),
    ("preview_svg", [["analyze", *P]], ["preview", *P, "--from-cache", "--frame", "5", "--output", "shot.svg"],
     "SVG of frame 5"),
    ("preview_pdf_clamped", [["analyze", *P]], ["preview", *P, "--from-cache", "--frame", "999", "--output",
                                                "out/last.pdf"], "frame beyond the end: the last one, PDF"),
    ("preview_negative", [["analyze", *P]], ["preview", *P, "--from-cache", "--frame", "-3"],
     "a negative frame: the first one"),
    ("preview_jpg", [["analyze", *P]], ["preview", *P, "--from-cache", "--output", "a.jpg"],
     "a format MDBridge does not offer (RN-11)"),
]


# --- the project file -------------------------------------------------------------------------------

def project_text(folder):
    p = demo_project("demo_dry.parm7", "demo_dry.nc", "results", 12, cleaning=False, mid_time_ps=120.0)
    p.analyses.mg.enabled = False
    p.analyses.stacking.enabled = False
    p.layout.preset = "large_plus_3"
    p.layout.panels = [PanelSpec("structure", [], "Structure"), PanelSpec("timeseries", ["rmsd/Complex"], "RMSD"),
                       PanelSpec("timeseries", ["hbonds/count"], ""),
                       PanelSpec("heatmap", ["contacts/per_residue_a"], "", {"max_ticks": 6})]
    p.output.timestamp_files = True
    p.output.pdb.enabled = False
    p.output.video.width, p.output.video.height, p.output.video.dpi, p.output.video.fps = 320, 180, 40, 12
    p.output.charts.dpi = 40
    p.output.charts.series = ["rmsd/Complex", "hbonds/count", "energy/total"]
    path = os.path.join(folder, PROJECT)
    p.save(path)
    with open(path, encoding="utf-8") as fh:
        text = fh.read()
    return SAVED_RE.sub('saved = "2026-11-25T00:00:00"', text)


# --- running ----------------------------------------------------------------------------------------

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
    return code, out.getvalue(), err.getvalue().replace(NETCDF_NOTICE, "")


def normalise(text, folder):
    for path in sorted({os.path.realpath(folder), os.path.abspath(folder)}, key=len, reverse=True):
        text = text.replace(path, "{work}")
    text = FINISHED_RE.sub("Finished in {seconds} s", text)
    return MANIFEST_RE.sub("{manifest}", CACHE_KEY_RE.sub("-{key}", STAMP_RE.sub("{stamp}", text)))


def mp4_info(path):
    probe = subprocess.run(["ffprobe", "-v", "error", "-count_frames", "-show_entries",
                            "stream=codec_type,codec_name,width,height,pix_fmt,r_frame_rate,avg_frame_rate,"
                            "time_base,nb_read_frames:format=format_name,nb_streams", "-of", "json", path],
                           capture_output=True, text=True, check=True)
    data = json.loads(probe.stdout)
    boxes = []
    with open(path, "rb") as fh:
        while True:
            head = fh.read(8)
            if len(head) < 8:
                break
            size, kind = struct.unpack(">I4s", head)
            if size == 1:
                size = struct.unpack(">Q", fh.read(8))[0]
                fh.seek(size - 16, 1)
            else:
                fh.seek(size - 8, 1)
            boxes.append(kind.decode("latin-1"))
            if size == 0:
                break
    streams = data.get("streams", [])
    return {"streams": streams, "format": data.get("format", {}), "boxes": boxes,
            "faststart": "moov" in boxes and "mdat" in boxes and boxes.index("moov") < boxes.index("mdat"),
            "audio_streams": sum(1 for st in streams if st.get("codec_type") == "audio")}


def file_info(full, rel, folder):
    with open(full, "rb") as fh:
        data = fh.read()
    info = {"magic": data[:8].hex(), "nonempty": len(data) > 0}
    if data[:8] == b"\x89PNG\r\n\x1a\n":
        info["size_px"] = list(struct.unpack(">II", data[16:24]))
    elif data[:5] == b"%PDF-":
        m = re.search(rb"/MediaBox \[ *([0-9. ]+)\]", data)
        info["size_pt"] = [float(v) for v in m.group(1).split()] if m else None
    elif b"<svg" in data[:2000]:
        m = re.search(rb'<svg[^>]*width="([0-9.]+)pt"[^>]*height="([0-9.]+)pt"', data)
        info["size_pt"] = [float(m.group(1)), float(m.group(2))] if m else None
    elif rel.endswith(".mp4"):
        info["video"] = mp4_info(full)
    elif rel.endswith("MDATK_params_log.toml"):
        doc = tomllib.loads(data.decode("utf-8"))
        run_ = doc.get("run", {})
        info["outputs"] = masked_log(run_.get("outputs"), folder)
        info["timestamp"] = masked_log(run_.get("timestamp"), folder)
    elif data[:4] == b"PK\x03\x04" or rel.endswith(".json"):
        info = {"cache": True}
    else:
        info["lines"] = data.count(b"\n")
    return info


def masked_log(entries, folder):
    def mask(v, key=""):
        if isinstance(v, dict):
            return {k: mask(x, k) for k, x in v.items()}
        if isinstance(v, list):
            return [mask(x) for x in v]
        if hasattr(v, "isoformat"):
            return "{time}"
        if isinstance(v, str):
            return "{time}" if ISO_RE.match(v) else normalise(v, folder)
        if key in ("duration_s",) and isinstance(v, (int, float)):
            return 0.0
        return v
    return mask(entries)


def run_case(name, before, argv, description, project):
    work = tempfile.mkdtemp(prefix="mdbridge-render-cli-")
    folder = os.path.join(work, "case")
    os.makedirs(folder)
    try:
        for dest, source in INPUTS.items():
            shutil.copyfile(source, os.path.join(folder, dest))
        with open(os.path.join(folder, PROJECT), "w", encoding="utf-8", newline="\n") as fh:
            fh.write(project)
        start = os.getcwd()
        os.chdir(folder)
        try:
            for b in before:
                code, _out, err = run(b)
                if code:
                    raise RuntimeError(f"{name}: {b} failed in MDATK: {err[-2000:]}")
            code, stdout, stderr = run(argv)
            files = {}
            for here, _dirs, names in os.walk(folder):
                for n in names:
                    full = os.path.join(here, n)
                    rel = os.path.relpath(full, folder).replace(os.sep, "/")
                    if rel in INPUTS or rel == PROJECT or (n.startswith(".") and n.endswith(".lock")):
                        continue
                    files[normalise(rel, folder)] = file_info(full, rel, folder)
            step_log = masked_log(read_log(os.path.join(folder, PROJECT)), folder)
        finally:
            os.chdir(start)
        rec = {"name": name, "description": description, "before": before, "argv": argv, "exit": code,
               "stdout": normalise(stdout, folder)}
        if TRACEBACK in stderr:
            rec["stderr"] = None
            rec["crash"] = normalise(stderr.strip().splitlines()[-1], folder)
        else:
            rec["stderr"] = normalise(stderr, folder)
        rec["files"] = dict(sorted(files.items()))
        rec["step_log"] = step_log
        return improvements_render.adjust("cli", name, rec)
    finally:
        shutil.rmtree(work, ignore_errors=True)


def main():
    if shutil.which("ffmpeg") is None or shutil.which("ffprobe") is None:
        sys.exit("make_render_cli_cases.py needs ffmpeg and ffprobe on PATH")
    os.makedirs(os.path.join(OUT, "cli"), exist_ok=True)
    tmp = tempfile.mkdtemp(prefix="mdbridge-render-cli-project-")
    try:
        project = project_text(tmp)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    cases = []
    for case in CASES:
        rec = run_case(*case, project)
        cases.append(rec)
        print(case[0], rec["exit"], rec.get("crash", ""), flush=True)
    rc.write_json(os.path.join(OUT, "cli", "render_cli_cases.json"),
                  {"improvements": improvements_render.applied(), "project_file": project, "cases": cases})


if __name__ == "__main__":
    main()
