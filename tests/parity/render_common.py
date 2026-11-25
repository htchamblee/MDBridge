"""Shared by the milestone-4 render generators (make_figure_cases.py, make_chart_cases.py,
make_molecule_cases.py, make_render_cli_cases.py; tests/data/README.md "Rendering").

Import it before numpy: it pins numpy's and OpenBLAS' kernels (as make_analysis_cases.py does), so
that what MDATK computes with them (view rotations, matplotlib's numerics, the software renderer's
float32 arithmetic, scipy.ndimage's Gaussian kernels) is the same on every x86-64 machine. Then
``setup(checkout)`` makes Matplotlib use only the fonts of its wheel (bundled_fonts_only), puts MDATK
on sys.path and patches it the way the analysis generators do:
MDBridge's independent tables (independent_tables.py), periodic geometry, analysis, selection and
CLI improvements, the residue-scoped ion rule on every Universe the pipeline loads, and MDBridge's
rendering improvements (improvements_render.py).

Also here:

* the input datasets of the figure and chart cases (``save_dataset`` / ``load_dataset``): every
  case is built from a dataset committed under tests/data/render/inputs/<name>/ (dataset.json and
  one <result>.npz per analysis result), so the C++ tests read exactly the data MDATK drew, without
  running the analyses; float64 arrays are rounded to 10 significant digits first, which keeps the
  committed bytes the same where a CPU's kernels move an aligned coordinate by an ulp;
* a deterministic .npz writer (members stored, not deflated: deflate output depends on the zlib
  build; fixed time stamps and permissions), readable by np.load and analysis::readNpz;
* ``dumps``: JSON with Python's float repr, non-finite floats as {"$float": "nan" / "inf" /
  "-inf"}; short objects and lists on one line, long objects one key per line (diffable, compact).
"""

from __future__ import annotations

import os
import platform
import sys

NATIVE = "--native" in sys.argv
if not NATIVE and "numpy" not in sys.modules:  # before numpy is imported (read when OpenBLAS / numpy load)
    os.environ["OPENBLAS_CORETYPE"] = "Nehalem"
    for _var in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ[_var] = "1"
    if platform.machine().lower() in ("x86_64", "amd64"):
        os.environ["NPY_ENABLE_CPU_FEATURES"] = "X86_V2"
os.environ["MPLBACKEND"] = "Agg"
os.environ.pop("DISPLAY", None)
os.environ.pop("WAYLAND_DISPLAY", None)

import io  # noqa: E402
import json  # noqa: E402
import logging  # noqa: E402
import math  # noqa: E402
import warnings  # noqa: E402
import zipfile  # noqa: E402

import numpy as np  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.normpath(os.path.join(HERE, "..", "data"))
DEMO_TOP = os.path.join(DATA, "topology", "demo_dry.parm7")
DEMO_NC = os.path.join(DATA, "io", "trajectories", "demo_dry.nc")
SYSTEMS = os.path.join(DATA, "analysis", "systems")

if HERE not in sys.path:
    sys.path.insert(0, HERE)


def args(usage: str) -> "tuple[str, str]":
    rest = [a for a in sys.argv[1:] if a != "--native"]
    if len(rest) != 2:
        sys.exit(f"usage: {usage}")
    return os.path.abspath(rest[0]), os.path.abspath(rest[1])


# --- MDATK, patched -------------------------------------------------------------------------------

LOG_LINES: list = []  # "LEVEL message" of every log record of MDATK's loggers since reset_log()


class _Collect(logging.Handler):
    def emit(self, record):
        LOG_LINES.append(f"{record.levelname} {record.getMessage()}")


def reset_log() -> None:
    LOG_LINES.clear()


_SETUP = []


def setup(checkout: str, collect: bool = True) -> None:
    """MDATK importable from `checkout`, patched like the analysis generators patch it. With
    `collect`, the records of MDATK's loggers go to LOG_LINES (and nowhere else); without, they
    propagate as in MDATK's command line (make_render_cli_cases.py)."""
    if _SETUP:
        return
    _SETUP.append(checkout)
    sys.path.insert(0, checkout)
    warnings.filterwarnings("ignore")
    logging.disable(logging.NOTSET)
    if collect:
        root = logging.getLogger("mdatk")
        root.setLevel(logging.INFO)
        root.propagate = False
        root.addHandler(_Collect())
    for noisy in ("MDAnalysis", "matplotlib", "PIL"):
        logging.getLogger(noisy).setLevel(logging.CRITICAL)

    import matplotlib

    matplotlib.use("Agg")
    matplotlib.rcdefaults()  # no matplotlibrc of the machine (MDATK uses the defaults)
    bundled_fonts_only()

    import independent_tables

    independent_tables.apply()  # MDAnalysis with MDBridge's independent tables (tests/data/README.md)
    import improvements_periodic

    improvements_periodic.apply_mdatk()
    import improvements_analyses

    improvements_analyses.apply_mdatk()
    import improvements_cli

    improvements_cli.apply()
    import improvements_readers_pdb  # noqa: F401  (MDBridge's reader warnings, as make_pipeline_cases.py)

    import mdatk.pipeline as pipeline_module
    from mdatk.io import loader

    original = loader.load_universe

    def load_universe(topology, trajectories=(), dt_ps=None):
        """load_universe with MDBridge's residue-scoped ion rule (as make_pipeline_cases.py loads)."""
        u = original(topology, trajectories, dt_ps)
        fmt = loader._TOPOLOGY_FORMATS.get(os.path.splitext(topology)[1].lower())
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            independent_tables.residue_rule_universe(u, topology, **({"topology_format": fmt} if fmt else {}))
        return u

    pipeline_module.load_universe = load_universe

    import improvements_render

    improvements_render.apply_mdatk()  # MDBridge's rendering improvements (RN-xx, docs/DIFFERENCES.md)
    import improvements_shading

    improvements_shading.apply_mdatk()  # residue shading, orientation, example data (RN-21 ... RN-25)
    import improvements_citations

    improvements_citations.apply_mdatk()  # MDBridge's citations (CIT-01, CIT-02; params and step logs)


def bundled_fonts_only() -> None:
    """Matplotlib finds only the fonts of its own wheel (DejaVu, pinned with it): a machine's own
    copies of DejaVu (other versions, other metrics) or other fonts never take part."""
    import matplotlib
    from matplotlib import font_manager as fm

    data = os.path.join(os.path.dirname(os.path.abspath(matplotlib.__file__)), "mpl-data")
    fm.fontManager.ttflist = [f for f in fm.fontManager.ttflist if os.path.abspath(f.fname).startswith(data)]
    fm.fontManager.afmlist = [f for f in fm.fontManager.afmlist if os.path.abspath(f.fname).startswith(data)]
    fm.fontManager._findfont_cached.cache_clear()
    fm._get_font.cache_clear()


# --- JSON -----------------------------------------------------------------------------------------

def plain(o):
    """JSON-able form: numpy scalars / arrays as Python values, tuples as lists, non-finite floats
    tagged."""
    if isinstance(o, dict):
        return {str(k): plain(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [plain(v) for v in o]
    if isinstance(o, np.ndarray):
        return plain(o.tolist())
    if isinstance(o, np.generic):
        return plain(o.item())
    if isinstance(o, float):
        if math.isnan(o):
            return {"$float": "nan"}
        if math.isinf(o):
            return {"$float": "inf" if o > 0 else "-inf"}
        return o
    return o


def unplain(o):
    if isinstance(o, dict):
        if set(o) == {"$float"}:
            return float(o["$float"])
        return {k: unplain(v) for k, v in o.items()}
    if isinstance(o, list):
        return [unplain(v) for v in o]
    return o


def _scalar(o) -> bool:
    return not isinstance(o, (dict, list)) or (isinstance(o, dict) and set(o) == {"$float"})


def _flat(o) -> bool:
    return isinstance(o, list) and all(_scalar(v) or (isinstance(v, list) and all(_scalar(w) for w in v)
                                                         and len(v) <= 4) for v in o)


ONE_LINE = 1000  # characters: shorter objects and lists are written on one line


def dumps(o, indent: int = 0) -> str:
    """Objects and lists of up to ONE_LINE characters on one line, longer objects one key per line,
    longer lists of scalars (or of short scalar lists: points, colours) on one line; floats as repr."""
    line = json.dumps(o, ensure_ascii=False, separators=(",", ":"), allow_nan=False)
    if len(line) <= ONE_LINE or _scalar(o) or _flat(o):
        return line
    pad = " " * indent
    if isinstance(o, dict):
        items = [f"{pad} {json.dumps(k, ensure_ascii=False)}: {dumps(v, indent + 1)}" for k, v in o.items()]
        return "{\n" + ",\n".join(items) + "\n" + pad + "}"
    return "[\n" + ",\n".join(f"{pad} {dumps(v, indent + 1)}" for v in o) + "\n" + pad + "]"


def write_json(path: str, obj) -> None:
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(dumps(plain(obj)) + "\n")


# --- .npz -----------------------------------------------------------------------------------------

def write_npz(path: str, arrays) -> None:
    """A numpy .npz archive (members stored, fixed time stamps): the same bytes on every machine."""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with zipfile.ZipFile(path, "w", zipfile.ZIP_STORED) as zf:
        for name, a in arrays.items():
            info = zipfile.ZipInfo(name + ".npy", date_time=(1980, 1, 1, 0, 0, 0))
            info.create_system = 3
            info.external_attr = 0o644 << 16
            buf = io.BytesIO()
            a = np.asarray(a)
            if a.dtype.byteorder == ">":
                a = a.astype(a.dtype.newbyteorder("<"))
            np.lib.format.write_array(buf, np.ascontiguousarray(a), allow_pickle=False)
            zf.writestr(info, buf.getvalue())


def read_npz(path: str) -> dict:
    with np.load(path, allow_pickle=False) as z:
        return {k: z[k] for k in z.files}


def round_sig(a: np.ndarray, digits: int = 10) -> np.ndarray:
    """float64 values rounded to `digits` significant digits (as float(f"{x:.{digits}g}")); other
    dtypes unchanged."""
    a = np.asarray(a)
    if a.dtype != np.float64 or a.size == 0:
        return a
    flat = [float(f"{x:.{digits}g}") for x in a.ravel().tolist()]
    return np.asarray(flat, dtype=np.float64).reshape(a.shape)


# --- input datasets -------------------------------------------------------------------------------

MASKS: list = []  # (path, placeholder) replaced before DATA -> {data} (add_mask)


def add_mask(path: str, placeholder: str) -> None:
    """Write `path` as `placeholder` (a generator's OUT_DIR as {data}/<its folder>)."""
    MASKS.append((os.path.abspath(path), placeholder))


def mask_text(s: str) -> str:
    for path, placeholder in MASKS:
        s = s.replace(path, placeholder)
    return s.replace(DATA, "{data}")


def _mask(o):
    if isinstance(o, dict):
        return {k: _mask(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_mask(v) for v in o]
    if isinstance(o, str):
        return mask_text(o)
    return o


def save_dataset(ds, folder: str) -> None:
    """tests/data/render/inputs/<name>/: dataset.json ({frames, time_ps, time_unit, results: [{name,
    key, params, meta, arrays (names in order)}]}) and <result>.npz per result (its arrays, float64
    ones rounded by round_sig). Cache keys are written as "{key}" (they hash input paths and
    mtimes)."""
    os.makedirs(folder, exist_ok=True)
    for name in os.listdir(folder):
        os.remove(os.path.join(folder, name))
    results = []
    for name, res in ds.results.items():
        results.append({"name": res.name, "key": "{key}", "params": _mask(plain(res.params)),
                        "meta": _mask(plain(res.meta)), "arrays": list(res.arrays)})
        write_npz(os.path.join(folder, f"{name}.npz"), {k: round_sig(v) for k, v in res.arrays.items()})
    write_json(os.path.join(folder, "dataset.json"),
               {"frames": ds.frames, "time_ps": round_sig(ds.time_ps), "time_unit": ds.time_unit,
                "results": results})


def load_dataset(folder: str):
    """The MDATK Dataset of save_dataset's files."""
    from mdatk.data import AnalysisResult, Dataset

    with open(os.path.join(folder, "dataset.json"), encoding="utf-8") as fh:
        doc = unplain(json.load(fh))
    ds = Dataset(np.asarray(doc["frames"], dtype=np.int64), np.asarray(doc["time_ps"], dtype=np.float64),
                 time_unit=doc["time_unit"])
    for r in doc["results"]:
        arrays = read_npz(os.path.join(folder, f"{r['name']}.npz"))
        ordered = {k: arrays[k] for k in r["arrays"]}
        meta = _unmask(r["meta"])
        ds.add(AnalysisResult(r["name"], r["key"], _unmask(r["params"]), ordered, meta))
    return ds


def _unmask(o):
    if isinstance(o, dict):
        return {k: _unmask(v) for k, v in o.items()}
    if isinstance(o, list):
        return [_unmask(v) for v in o]
    if isinstance(o, str):
        return o.replace("{data}", DATA)
    return o


def project_dict(project) -> dict:
    """Project.to_dict() as tagged JSON, paths of tests/data as {data}."""
    return _mask(plain(project.to_dict()))


def project_from(d: dict):
    from mdatk.config import Project

    return Project.from_dict(_unmask(unplain(d)), None)
