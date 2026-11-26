"""What Matplotlib 3.11.2 (as MDATK v0.0.4.1 uses it) does with the IBM Plex fonts MDBridge bundles
(data/fonts) when it is configured as MDBridge draws (improvements_fonts.py: the IBM Plex Serif and
Plex Mono faces only, Plex Serif first in the sans-serif / serif lists, Plex Mono first in the
monospace list, IBM Plex Math as the fallback font), for
src/render/core/font.cpp and text.cpp (tests/render/test_text.cpp).

    /opt/mdbridge-venv/bin/python -I -B make_text_cases.py OUT_DIR

Writes OUT_DIR/text_cases.json with the sections:

  fonts      every bundled file (the six figure faces: four Plex Serif, two Plex Mono for
             sequences and the interface's fixed-width text; the fallback) with its role: file
             name, sha256 (IBM's release file, checked by improvements_fonts.path), family and
             style as FreeType reports them, the name table's version, units per em, ascender,
             descender, glyph count. Matplotlib's FreeType version is recorded too.
  choice     the face Matplotlib's font manager picks (font_manager.findfont) among the six
             figure faces (in data/fonts' registry order; the machine's fonts left out, so the
             result is the same everywhere) for font properties: family lists (generic names,
             the Plex and DejaVu names, unknown families), weight (names and numbers) and style
             (normal, italic, oblique).
  metrics    RendererAgg.get_text_width_height_descent (the extents Text layout uses) for many
             strings (MDATK's labels, tick labels with Unicode minus, superscripts, ion charges
             from the fallback font, Angstrom, Greek, ligatures, kerning pairs, empty and blank
             strings) on every figure face at several font sizes and dpi, with text.hinting
             "no_hinting" (FreeType's unhinted outlines, which MDBridge reproduces exactly: width /
             height / descent and the glyph x positions in 26.6 fixed point from FT2Font.set_text
             with the fallback font) and with Matplotlib's default hinting (what MDATK's figures
             measure; MDBridge differs within a bound, docs/DIFFERENCES.md).

Matplotlib's rcParams are its defaults (no matplotlibrc of the machine) with MDBridge's font lists.
The output depends only on Matplotlib's wheel (its FreeType and HarfBuzz builds) and the fonts: the
same on every machine.
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

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import improvements_fonts  # noqa: E402

improvements_fonts.apply()

from matplotlib import font_manager as fm  # noqa: E402
from matplotlib import ft2font  # noqa: E402
from matplotlib.backends.backend_agg import RendererAgg  # noqa: E402

OUT = os.path.abspath(sys.argv[1])
FILES = [f for f, *_ in improvements_fonts.FACES]  # the figure faces, registry order
ROLES = ([(f, "sequences, interface fixed-width" if fam == improvements_fonts.SEQUENCE_FAMILY else "figures, interface")
          for f, fam, *_ in improvements_fonts.FACES] + [(improvements_fonts.FALLBACK[0], "fallback")])

STRINGS = [
    "", " ", "  ", "lp", "Hello World", "Time (ns)", "Time (µs)", "RMSD (Å)", "RMSD to frame 0 (Å)",
    "frame 120", "Radius of gyration (nm)", "fraction of frames", "count", "Distance to Mg²⁺ (Å)",
    "Mg²⁺", "K⁺", "Na⁺", "Ions", "Chain A", "Chain 12", "A:12", "B:1043", "No data: no frames",
    "−0.5", "−12.75", "1e−3", "×10³", "0", "0123456789", "100", "2.5", "Δt = 2 ps", "α β γ",
    "AUGC", "G", "U", "WAVE Tj", "To AV", "Ty LT", "fi fl ffi office", "jjj", "yjg", "Water",
    "TeV", "LAVA", "Energy (kJ/mol)", "H-bonds", "contacts", "Base pairs", "MDATK v0.0.4.1",
    "€ ‰ ±", "…", "[1, 2]", "{x}", "50 %", "#", "_",
    # characters IBM Plex Serif lacks, from the fallback font (IBM Plex Math): MDBridge's ion legend
    # entries (RN-08), fallback runs between Plex Serif runs, a ligature and a kerning pair next to
    # them, Plex Math alone
    "Cl⁻", "Ca²⁺", "Zn²⁺ and Mn²⁺", "Na⁺ / K⁺", "AV⁺fi", "⁺⁻", "⁺", "x⁺⁻y (Å)", "Ions: Mg²⁺, Cl⁻",
    "I l 1 | Il1",
    # sequence text (IBM Plex Mono faces): one-letter sequences, residue labels, dot-bracket
    "GCAUGCAUGC", "A:G15", "((((...))))..", "ACGU acgu",
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
    for name, role in ROLES:
        path = improvements_fonts.path(name)
        font = ft2font.FT2Font(path)
        sfnt = font.get_sfnt()
        version = sfnt.get((3, 1, 0x409, 5)) or sfnt.get((1, 0, 0, 5))
        version = version.decode("utf-16-be" if (3, 1, 0x409, 5) in sfnt else "latin-1")
        out.append({"file": name, "role": role, "sha256": sha256(path), "family": font.family_name,
                    "style": font.style_name, "version": version, "units_per_em": font.units_per_EM,
                    "ascender": font.ascender, "descender": font.descender, "glyphs": font.num_glyphs})
    return out


def choice_section():
    manager = fm.fontManager  # the figure faces only, in registry order (improvements_fonts.apply)
    assert [os.path.basename(e.fname) for e in manager.ttflist] == FILES
    families = [["sans-serif"], ["IBM Plex Serif"], ["monospace"], ["IBM Plex Mono"], ["DejaVu Sans"],
                ["Unknown Family"], ["DejaVu Sans Mono", "serif"],
                ["sans"], ["serif"], ["Unknown Family", "monospace"], ["Arial", "IBM Plex Serif"],
                ["ibm plex serif", "sans-serif"], ["cursive"], ["fantasy"], ["IBM Plex Math"]]
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
        path = improvements_fonts.path(name)
        font = fm.get_font([path, improvements_fonts.fallback_path()])  # as the renderer gets it
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
