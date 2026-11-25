# Covalent radii

- **`cordero2008.tsv`**: Z, element symbol and single-bond covalent radius (Å) for Z 1-96, from
  B. Cordero, V. Gómez, A. E. Platero-Prats, M. Revés, J. Echeverría, E. Cremades, F. Barragán and
  S. Alvarez, *Covalent radii revisited*, Dalton Trans. (2008) 2832-2838, doi:10.1039/b801115j,
  Table 2. Carbon is the sp3 value; Mn, Fe and Co are the low-spin values (the other values are
  noted). Used by the distance bond guesser of the analyses (AN-14, docs/IMPROVEMENTS.md):
  `tools/gen_covalent_radii.py` writes `src/analysis/covalent_radii.inc` from this file.

**Licence:** radii are published scientific data (facts); no copyrightable expression of the
article is reproduced. Independent of MDAnalysis (whose bond guesser uses van der Waals radii).
