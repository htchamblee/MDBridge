"""Draw MDBridge's application icon: two atoms (the residue-group colours) joined by a bridge.

    python packaging/make_icon.py [OUT_DIR]      (default: data/icons)

Writes mdbridge.png (256 px) and mdbridge.ico (16-256 px). Both are committed, so a build does not
need this script or Pillow; run it again only to change the icon.
"""

from __future__ import annotations

import os
import sys

from PIL import Image, ImageDraw, ImageFilter

SIZE = 1024  # drawn large, then scaled down (anti-aliasing)
BACKGROUND = (26, 30, 39)
BLUE = (42, 120, 214)
ORANGE = (235, 104, 52)
BRIDGE = (226, 232, 240)


def sphere(draw: ImageDraw.ImageDraw, centre, radius, colour) -> None:
    x, y = centre
    draw.ellipse([x - radius, y - radius, x + radius, y + radius], fill=colour)
    for k in range(12):  # soft highlight towards the upper left
        r = radius * (0.55 - k * 0.04)
        shade = tuple(min(255, int(c + (255 - c) * (0.08 + k * 0.05))) for c in colour)
        ox, oy = x - radius * 0.28, y - radius * 0.28
        draw.ellipse([ox - r, oy - r, ox + r, oy + r], fill=shade)


def icon() -> Image.Image:
    img = Image.new("RGBA", (SIZE, SIZE), (0, 0, 0, 0))
    draw = ImageDraw.Draw(img)
    draw.rounded_rectangle([40, 40, SIZE - 40, SIZE - 40], radius=200, fill=BACKGROUND)
    left, right, base = (300, 640), (724, 640), 640
    # the bridge: an arc from atom to atom, with hangers down to the deck
    draw.arc([left[0], 300, right[0], 980], start=200, end=340, fill=BRIDGE, width=44)
    draw.line([left, right], fill=BRIDGE, width=36)
    for x in (412, 512, 612):
        top = 360 if x == 512 else 395
        draw.line([(x, top), (x, base)], fill=BRIDGE, width=22)
    sphere(draw, left, 150, BLUE)
    sphere(draw, right, 150, ORANGE)
    return img.filter(ImageFilter.SMOOTH)


def main(argv=None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    out = argv[0] if argv else os.path.join(root, "data", "icons")
    os.makedirs(out, exist_ok=True)
    big = icon()
    big.resize((256, 256), Image.LANCZOS).save(os.path.join(out, "mdbridge.png"), optimize=True)
    big.resize((256, 256), Image.LANCZOS).save(
        os.path.join(out, "mdbridge.ico"), sizes=[(s, s) for s in (16, 24, 32, 48, 64, 128, 256)])
    print(f"wrote {out}/mdbridge.png and mdbridge.ico")
    return 0


if __name__ == "__main__":
    sys.exit(main())
