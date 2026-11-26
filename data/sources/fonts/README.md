# The bundled fonts (data/fonts)

MDBridge draws its interface and every figure, chart and video text in fonts it carries itself,
compiled into the programs as Qt resources (`src/render/CMakeLists.txt`). The files in `data/fonts`
are unmodified TrueType files of IBM's official IBM Plex releases (https://github.com/IBM/plex),
licensed under the SIL Open Font License 1.1 (`packaging/licenses/IBM-Plex-OFL-1.1.txt`, the
releases' `LICENSE.txt`: Copyright © 2017 IBM Corp. with Reserved Font Name "Plex"). Because the
files are not modified, the reserved font name needs no change; a modified copy would have to be
renamed (OFL condition 3).

| File | Role | Release tag | Release asset | Path in the asset | sha256 of the file |
|---|---|---|---|---|---|
| `IBMPlexSerif-Regular.ttf` | figures, interface | `@ibm/plex-serif@2.0.0` | `ibm-plex-serif.zip` | `ibm-plex-serif/fonts/complete/ttf/` | `6c490ccd38e856c674e0d163f3f40742a75058b7b935998dbfbc1745ba54e474` |
| `IBMPlexSerif-Bold.ttf` | figures, interface | `@ibm/plex-serif@2.0.0` | `ibm-plex-serif.zip` | same | `5fa66833dc80adb6252c793adb555716ec10833d2b1e5b0235b68c72cd69adf0` |
| `IBMPlexSerif-Italic.ttf` | figures, interface | `@ibm/plex-serif@2.0.0` | `ibm-plex-serif.zip` | same | `d28aa63cf618439712d6b459c5649af1d6d789cb911994f2357f1ce44128f163` |
| `IBMPlexSerif-BoldItalic.ttf` | figures, interface | `@ibm/plex-serif@2.0.0` | `ibm-plex-serif.zip` | same | `e33c9d0a692bb26c9c1ece4a35fd2a351317f5a60009d0e9945ae30d0297bf42` |
| `IBMPlexMath-Regular.ttf` | fallback for missing glyphs (figures and interface) | `@ibm/plex-math@1.1.0` | `ibm-plex-math.zip` | `ibm-plex-math/fonts/complete/ttf/` | `10d07a23cb322725f83ca51748a4c7a61cf14caf2c4d880d793c65d3a1402533` |
| `IBMPlexMono-Regular.ttf` | interface fixed-width text (logs, multi-line fields) | `@ibm/plex-mono@2.5.0` | `ibm-plex-mono.zip` | `ibm-plex-mono/fonts/complete/ttf/` | `7c6fbddca4b700be918f5f6183d9bd4464fa427fe435f0b480d77fe2bb8c5a43` |
| `IBMPlexMono-Bold.ttf` | interface fixed-width text | `@ibm/plex-mono@2.5.0` | `ibm-plex-mono.zip` | same | `74e5eedcfa4596497d34e19023cabdabd3a8c852b903007a5654a59591a72ffb` |

The release assets (downloaded 2026-11-26):

| Asset | URL | sha256 |
|---|---|---|
| `ibm-plex-serif.zip` | https://github.com/IBM/plex/releases/download/%40ibm%2Fplex-serif%402.0.0/ibm-plex-serif.zip | `c006b095c47b919c6c8d78319e9e628745033f5db94507c1c9dfe4b35dd51f1e` |
| `ibm-plex-math.zip` | https://github.com/IBM/plex/releases/download/%40ibm%2Fplex-math%401.1.0/ibm-plex-math.zip | `d85ed404394ced3a79a519af24611acdee9cc0483363c07cd5ac0656c06db42a` |
| `ibm-plex-mono.zip` | https://github.com/IBM/plex/releases/download/%40ibm%2Fplex-mono%402.5.0/ibm-plex-mono.zip | `6d23f01257663d8cc49a0d64c22ced630b79e0e2a0ac08a0da86e9a38bbc481c` |

The same checksums of the files are pinned in `src/render/core/text.cpp` (the programs refuse other
files) and `tests/parity/improvements_fonts.py` (the parity generators give Matplotlib the same
files).

## Why these faces

- **IBM Plex Serif** (Regular, Bold, Italic, Bold Italic): the owner's choice for all text, in the
  interface and in figures, because its I, l and 1 cannot be mistaken for one another (DejaVu Sans,
  Matplotlib's default and MDATK's figure font, draws I and l as plain strokes). Figures use
  Regular, Bold (the sequence panel's letters) and Italic (chain notes); Bold Italic is what a bold
  italic text gets (Matplotlib's font choice, `chooseFace`) and what Qt uses for bold italic interface text.
- **IBM Plex Math** (Regular, the only face): the one fallback for the characters Plex Serif has
  no glyph for. Of the characters MDBridge's figures and interface can show, these are SUPERSCRIPT
  PLUS SIGN (U+207A) and SUPERSCRIPT MINUS (U+207B), in the ion labels (Mg²⁺, Na⁺, Cl⁻ ...);
  Plex Math has both, and every other character that MDBridge's figures and interface can show
  is in Plex Serif itself (checked over all non-ASCII characters in `src/`, the text cases'
  strings and MDATK's figure labels; `tests/render/test_text.cpp`). No character was found that
  neither font has.
- **IBM Plex Mono** (Regular, Bold): fixed-width text needs a monospace font (the log view's
  aligned columns, TOML and command fields); Plex Mono has the same unambiguous I / l / 1.
  Figures do not use it: MDATK's "monospace" sequence letters are single centred letters, drawn in
  Plex Serif like the rest (docs/DIFFERENCES.md FONT-01). Plex Mono lacks Greek (ε, σ appear in
  some labels, never in fixed-width text); Plex Math is its fallback too.

DejaVu (the files Matplotlib 3.11.2 ships) was bundled until this change and is gone: nothing needs
it any more.
