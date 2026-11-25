"""MDBridge's colour improvements (docs/DIFFERENCES.md "Colours and themes"), for
make_color_cases.py: the expected values where MDBridge deliberately differs from MDATK /
Matplotlib.

  RN-18  theme.ordinal(n) for more than 5 classes: n - 1 distinct steps along the 4-colour ramp,
         the lookup table of LinearSegmentedColormap.from_list(ramp, N = n - 1) (MDATK's
         ListedColormap kept the 4 ramp colours, so with vmin -0.5 / vmax n - 0.5 several classes
         shared a colour)
  RN-12  an RGB(A) tuple with a NaN component is not a colour ("RGBA values should be within 0-1
         range"); Matplotlib's range test lets NaN through
  RN-13  an unknown or unbundled colormap name: MDBridge's message listing the bundled maps (owner
         decision D6), with a suggestion for a name that differs only in case
"""

from __future__ import annotations

import math

from matplotlib.colors import LinearSegmentedColormap, to_rgba

# owner decision D6 (tools/gen_mpl_tables.py BUNDLED_MAPS, the same order)
BUNDLED_MAPS = ["tab10", "tab20", "tab20b", "tab20c", "Set1", "Set2", "Set3", "Pastel1", "Pastel2", "Dark2",
                "Accent", "Paired", "viridis", "plasma", "inferno", "magma", "cividis"]

ORDINAL_RAMP = {"light": ["#86b6ef", "#3987e5", "#1c5cab", "#0d366b"],
                "dark": ["#184f95", "#3987e5", "#86b6ef", "#cde2fb"]}
ORDINAL_ZERO = {"light": "#e1e0d9", "dark": "#2c2c2a"}


def ordinal_colors(theme_name: str, n: int, mdatk_colors):
    """The colours of MDBridge's ordinal(n) (MDATK's for n <= 5)."""
    ramp = ORDINAL_RAMP["dark" if theme_name == "dark" else "light"]
    if n - 1 <= len(ramp):
        return mdatk_colors
    zero = ORDINAL_ZERO["dark" if theme_name == "dark" else "light"]
    along = LinearSegmentedColormap.from_list("ordinal", ramp, N=n - 1)
    along._init()
    return [list(to_rgba(zero))] + [[float(v) for v in row] for row in along._lut[: n - 1]]


def tuple_result(components, alpha, matplotlib_result):
    """Matplotlib's to_rgba of a tuple, but NaN components are an error."""
    if isinstance(matplotlib_result, list) and any(math.isnan(v) for v in matplotlib_result):
        return {"error": "RGBA values should be within 0-1 range"}
    return matplotlib_result


def _repr(text: str) -> str:
    return repr(text)


def colormap_problem(name: str):
    """MDBridge's message for a colormap name that is not bundled, or None."""
    base = name[:-2] if len(name) > 2 and name.endswith("_r") else name
    if base in BUNDLED_MAPS:
        return None
    message = (f"{_repr(name)} is not a bundled colormap; choose from {', '.join(BUNDLED_MAPS)} "
               '(or one of them reversed, with "_r")')

    def lower(text):  # ASCII letters only, as the C++ code
        return "".join(chr(ord(c) + 32) if "A" <= c <= "Z" else c for c in text)

    for candidate in BUNDLED_MAPS:
        for spelled in (candidate, candidate + "_r"):
            if lower(spelled) == lower(name):
                return message + f"; did you mean {_repr(spelled)}?"
    return message
