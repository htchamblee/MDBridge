"""MDBridge's figure fonts in Matplotlib (docs/DIFFERENCES.md FONT-01, FONT-02), in the style of the
other improvements_*.py modules: MDATK draws its text in the DejaVu Sans of Matplotlib's wheel;
MDBridge draws every figure, chart and video text in the IBM Plex Serif it bundles (data/fonts),
with IBM Plex Math as the one fallback for the characters Plex Serif lacks. ``apply()`` makes
Matplotlib do the same, so that the expected data of the render generators (text metrics, layouts,
tick label sizes) is what MDBridge draws:

* Matplotlib's font list is the four IBM Plex Serif faces of data/fonts only (Regular, Bold,
  Italic, Bold Italic, in that order: MDBridge's registry, src/render/core/text.cpp), checked
  against the sha256 of IBM's release (data/sources/fonts/README.md); the machine's fonts and the
  wheel's DejaVu files never take part.
* IBM Plex Serif comes first in font.sans-serif, font.serif and font.monospace (rcParams and
  rcParamsDefault, so the "default" chart style keeps it) and is the font manager's default family:
  every family a text asks for (MDATK's default sans-serif, its monospace sequence letters, unknown
  names) is drawn in Plex Serif, the face chosen by Matplotlib's style and weight scoring.
* Every font Matplotlib's Agg renderer, text-to-path and PDF backends look up for font properties
  (FontManager._find_fonts_by_props) gets data/fonts/IBMPlexMath-Regular.ttf as its next font, so
  FT2Font takes the characters Plex Serif has no glyph for (superscript + and -, ...) from Plex Math,
  per cluster, as it falls back between fonts (ft2font.cpp FT2Font::layout); MDBridge's
  Font::layout does the same (src/render/core/font.cpp).

Standard library only at import; Matplotlib is imported by ``apply()``.
"""

from __future__ import annotations

import hashlib
import os

HERE = os.path.dirname(os.path.abspath(__file__))
FONTS = os.path.normpath(os.path.join(HERE, "..", "..", "data", "fonts"))

# (file, family, style, sha256 of IBM's release file): the figure faces in registry order
FACES = [
    ("IBMPlexSerif-Regular.ttf", "IBM Plex Serif", "Regular",
     "6c490ccd38e856c674e0d163f3f40742a75058b7b935998dbfbc1745ba54e474"),
    ("IBMPlexSerif-Bold.ttf", "IBM Plex Serif", "Bold",
     "5fa66833dc80adb6252c793adb555716ec10833d2b1e5b0235b68c72cd69adf0"),
    ("IBMPlexSerif-Italic.ttf", "IBM Plex Serif", "Italic",
     "d28aa63cf618439712d6b459c5649af1d6d789cb911994f2357f1ce44128f163"),
    ("IBMPlexSerif-BoldItalic.ttf", "IBM Plex Serif", "Bold Italic",
     "e33c9d0a692bb26c9c1ece4a35fd2a351317f5a60009d0e9945ae30d0297bf42"),
]
FALLBACK = ("IBMPlexMath-Regular.ttf", "IBM Plex Math", "Regular",
            "10d07a23cb322725f83ca51748a4c7a61cf14caf2c4d880d793c65d3a1402533")
# the interface's fixed-width faces (not used in figures; checked by make_text_cases.py)
INTERFACE_MONO = [
    ("IBMPlexMono-Regular.ttf", "IBM Plex Mono", "Regular",
     "7c6fbddca4b700be918f5f6183d9bd4464fa427fe435f0b480d77fe2bb8c5a43"),
    ("IBMPlexMono-Bold.ttf", "IBM Plex Mono", "Bold",
     "74e5eedcfa4596497d34e19023cabdabd3a8c852b903007a5654a59591a72ffb"),
]
FAMILY = "IBM Plex Serif"
GENERIC = ("font.sans-serif", "font.serif", "font.monospace")

_APPLIED: list = []


def path(file: str) -> str:
    """A bundled font file, checked against IBM's release."""
    p = os.path.join(FONTS, file)
    expected = {f: h for f, _fam, _style, h in FACES + [FALLBACK] + INTERFACE_MONO}[file]
    with open(p, "rb") as fh:
        actual = hashlib.sha256(fh.read()).hexdigest()
    if actual != expected:
        raise RuntimeError(f"{p} is not IBM's release file (sha256 {actual})")
    return p


def fallback_path() -> str:
    return path(FALLBACK[0])


def apply() -> None:
    """Matplotlib with MDBridge's figure fonts (idempotent). Call after matplotlib.rcdefaults()."""
    import matplotlib as mpl
    from matplotlib import font_manager as fm

    if not _APPLIED:
        _APPLIED.append(True)
        faces = [path(f) for f, *_ in FACES]
        math = fallback_path()
        manager = fm.fontManager
        manager.ttflist = []
        manager.afmlist = []
        for p in faces:
            manager.addfont(p)
        # addfont lists a file once per family name it advertises; only the figure faces take part
        manager.ttflist = [e for e in manager.ttflist if os.path.abspath(e.fname) in faces]
        manager.defaultFamily = {"ttf": FAMILY, "afm": FAMILY}
        original = fm.FontManager._find_fonts_by_props

        def find_fonts_by_props(self, prop, *args, **kwargs):
            """Matplotlib's font list for `prop`, with IBM Plex Math as the next font (FONT-02)."""
            found = list(original(self, prop, *args, **kwargs))
            if not any(os.path.abspath(str(getattr(f, "path", f))) == math for f in found):
                found.append(math)
            return found

        fm.FontManager._find_fonts_by_props = find_fonts_by_props
    for rc in (mpl.rcParams, mpl.rcParamsDefault):
        for key in GENERIC:
            names = list(rc[key])
            if not names or names[0] != FAMILY:
                rc[key] = [FAMILY] + [n for n in names if n != FAMILY]
    fm.fontManager._findfont_cached.cache_clear()
    fm._get_font.cache_clear()
