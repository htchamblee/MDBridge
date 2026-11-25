"""What Matplotlib 3.11.2 (as MDATK v0.0.4.1 uses it) does with the DejaVu fonts MDBridge bundles
(data/fonts), for src/render/core/font.cpp and text.cpp (tests/render/test_text.cpp).

    /opt/mdbridge-venv/bin/python -I -B make_text_cases.py OUT_DIR

Writes OUT_DIR/text_cases.json with the sections:

  fonts      Matplotlib's own copies of the five bundled faces (mpl-data/fonts/ttf): file name,
             sha256 (the committed data/fonts files must have the same), family and style as
             FreeType reports them, the name table's version, units per em, ascender, descender,
             glyph count. Matplotlib's FreeType version is recorded too.
  choice     the face Matplotlib's font manager picks (font_manager.findfont) among the five
             bundled faces (in data/fonts' order; the machine's fonts left out, so the result is the
             same everywhere) for font properties: family lists (generic names, the DejaVu names,
             unknown families), weight (names and numbers) and style (normal, italic, oblique).
  metrics    RendererAgg.get_text_width_height_descent (the extents Text layout uses) for many
             strings (MDATK's labels, tick labels with Unicode minus, superscripts, Angstrom, Greek,
             ligatures, kerning pairs, empty and blank strings) on every face at several font sizes
             and dpi, with text.hinting "no_hinting" (FreeType's unhinted outlines, which MDBridge
             reproduces exactly: width / height / descent and the glyph x positions in 26.6 fixed
             point from FT2Font.set_text) and with Matplotlib's default hinting (what MDATK's
             figures measure; MDBridge differs within a bound, docs/DIFFERENCES.md).

Matplotlib's rcParams are its defaults (no matplotlibrc of the machine). The output depends only on
Matplotlib's wheel (its FreeType and HarfBuzz builds) and the fonts: the same on every machine.
"""

import hashlib
import json
import logging
import os
import sys
import tempfile

os.environ["MPLBACKEND"] = "Agg"
os.environ["MPLCONFIGDIR"] = tempfile.mkdtemp(prefix="mdbridge-mpl-")  # no user matplotlibrc

if len(sys.argv) != 2:
    sys.exit("usage: make_text_cases.py OUT_DIR")

import matplotlib as mpl  # noqa: E402

mpl.use("Agg")
mpl.rcdefaults()
logging.getLogger("matplotlib").setLevel(logging.ERROR)  # findfont's fallback notes

from matplotlib import font_manager as fm  # noqa: E402
from matplotlib import ft2font  # noqa: E402
from matplotlib.backends.backend_agg import RendererAgg  # noqa: E402

OUT = os.path.abspath(sys.argv[1])
TTF = os.path.join(mpl.get_data_path(), "fonts", "ttf")
FILES = ["DejaVuSans.ttf", "DejaVuSans-Bold.ttf", "DejaVuSans-Oblique.ttf", "DejaVuSansMono.ttf",
         "DejaVuSansMono-Bold.ttf"]

STRINGS = [
    "", " ", "  ", "lp", "Hello World", "Time (ns)", "Time (µs)", "RMSD (Å)", "RMSD to frame 0 (Å)",
    "frame 120", "Radius of gyration (nm)", "fraction of frames", "count", "Distance to Mg²⁺ (Å)",
    "Mg²⁺", "K⁺", "Na⁺", "Ions", "Chain A", "Chain 12", "A:12", "B:1043", "No data: no frames",
    "−0.5", "−12.75", "1e−3", "×10³", "0", "0123456789", "100", "2.5", "Δt = 2 ps", "α β γ",
    "AUGC", "G", "U", "WAVE Tj", "To AV", "Ty LT", "fi fl ffi office", "jjj", "yjg", "Water",
    "TeV", "LAVA", "Energy (kJ/mol)", "H-bonds", "contacts", "Base pairs", "MDATK v0.0.4.1",
    "€ ‰ ±", "…", "[1, 2]", "{x}", "50 %", "#", "_",
]
# (points, dpi): MDATK's font sizes (10 by default, x 0.6 ... 1.3) at the dpi of videos, charts and
# previews, and two non-integer dpi (FreeType takes the integer part)
SIZES = [(10.0, 100.0), (11.0, 150.0), (8.0, 72.0), (6.0, 200.0), (13.0, 100.0), (7.0, 96.0),
         (8.5, 120.0), (9.0, 300.0), (10.0, 72.7), (12.0, 150.5)]


def sha256(path):
    with open(path, "rb") as fh:
        return hashlib.sha256(fh.read()).hexdigest()


def fonts_section():
    out = []
    for name in FILES:
        path = os.path.join(TTF, name)
        font = ft2font.FT2Font(path)
        sfnt = font.get_sfnt()
        version = sfnt.get((3, 1, 0x409, 5)) or sfnt.get((1, 0, 0, 5))
        version = version.decode("utf-16-be" if (3, 1, 0x409, 5) in sfnt else "latin-1")
        out.append({"file": name, "sha256": sha256(path), "family": font.family_name, "style": font.style_name,
                    "version": version, "units_per_em": font.units_per_EM, "ascender": font.ascender,
                    "descender": font.descender, "glyphs": font.num_glyphs})
    return out


def choice_section():
    manager = fm.FontManager()
    # the bundled faces only, in data/fonts' order (MDBridge's registry): what Matplotlib's scoring
    # picks among them (the machine's fonts and Matplotlib's other DejaVu faces left out)
    entries = {os.path.basename(e.fname): e for e in manager.ttflist
               if os.path.dirname(os.path.abspath(e.fname)) == os.path.abspath(TTF)}
    manager.ttflist = [entries[name] for name in FILES]
    manager.afmlist = []
    families = [["sans-serif"], ["DejaVu Sans"], ["monospace"], ["DejaVu Sans Mono"], ["Unknown Family"],
                ["sans"], ["serif"], ["Unknown Family", "monospace"], ["Arial", "DejaVu Sans Mono"],
                ["dejavu sans mono", "sans-serif"], ["cursive"]]
    weights = ["normal", "bold", "light", "ultralight", "medium", "semibold", "demibold", "demi", "heavy",
               "black", "regular", "book", "roman", "extra bold", 100, 300, 400, 500, 549, 550, 551, 600, 700,
               900]
    out = []
    for family in families:
        for weight in weights:
            for style in ("normal", "italic", "oblique"):
                prop = fm.FontProperties(family=family, weight=weight, style=style)
                path = manager.findfont(prop, fallback_to_default=True, rebuild_if_missing=False)
                out.append({"family": family, "weight": weight, "style": style, "file": os.path.basename(path)})
    return out


def metrics_section():
    out = []
    k = 0
    for name in FILES:
        path = os.path.join(TTF, name)
        font = ft2font.FT2Font(path)
        for s in STRINGS:
            for j in range(4):  # four of the sizes per string, cycling through all of them
                pt, dpi = SIZES[(k + 3 * j) % len(SIZES)]
                renderer = RendererAgg(10, 10, dpi)
                prop = fm.FontProperties(fname=path, size=pt)
                with mpl.rc_context({"text.hinting": "no_hinting"}):
                    w, h, d = renderer.get_text_width_height_descent(s, prop, ismath=False)
                hw, hh, hd = renderer.get_text_width_height_descent(s, prop, ismath=False)
                font.clear()
                font.set_size(pt, dpi)
                xy = font.set_text(s, 0.0, flags=ft2font.LoadFlags.NO_HINTING)
                xs = [int(v) for v in xy.reshape(-1, 2)[:, 0]] if len(xy) else []
                out.append({"file": name, "text": s, "points": pt, "dpi": dpi,
                            "width": int(round(w * 64)), "height": int(round(h * 64)), "descent": int(round(d * 64)),
                            "x": xs, "hinted": [hw, hh, hd]})
            k += 1
    return out


def main():
    data = {
        "matplotlib": mpl.__version__,
        "freetype": ft2font.__freetype_version__,
        "fonts": fonts_section(),
        "choice": choice_section(),
        "metrics": metrics_section(),
    }
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, "text_cases.json"), "w", encoding="utf-8", newline="\n") as fh:
        fh.write("{\n")
        keys = list(data)
        for i, key in enumerate(keys):
            value = data[key]
            sep = "," if i + 1 < len(keys) else ""
            if isinstance(value, list):
                fh.write(f' "{key}": [\n')
                fh.write(",\n".join("  " + json.dumps(v, ensure_ascii=False, separators=(",", ":")) for v in value))
                fh.write(f"\n ]{sep}\n")
            else:
                fh.write(f' "{key}": {json.dumps(value, ensure_ascii=False)}{sep}\n')
        fh.write("}\n")


if __name__ == "__main__":
    main()
