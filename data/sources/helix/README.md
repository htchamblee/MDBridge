# Helix geometry: the Arnott A-RNA fibre monomer

| File | Origin | Retrieved | Licence |
|---|---|---|---|
| `arnott_arna.tsv` | numbers read from the `$ARNA` block ("ARNOTT ARNA MONOMER UNIT", lines 125-184) of AmberTools' `dat/leap/parm/nucgen.dat`: one row per line of the block with its line number, the file's atom name, the PDB v3 name, r (Angstrom), phi (degrees) and z (Angstrom) as the file writes them | 2026-11-26 | facts (published coordinates with their provenance; no file copied) |

The coordinates are the A-RNA fibre model of Arnott, Hukins & Dover, *Optimised parameters for RNA
double-helices*, Biochem. Biophys. Res. Commun. 48 (1972) 1392-1399, doi:10.1016/0006-291X(72)90867-4
(citation key `arnott1972`), the monomer of NAB's and LEaP's nucgen; Chandrasekaran & Arnott tabulate
the fibre models in Landolt-Börnstein, New Series VII/1b (1989). nucgen.dat ships with AmberTools,
whose components are under several licences (GPL, LGPL and others); MDBridge does not rely on that
licence: as for parm10.dat (`../forcefields/README.md`), only the numbers of the one block are kept,
with file and line, and no other content (other blocks, comments, layout) is copied. MDATK 0.0.4.1
embeds the same numbers (`builder.py` `_BACKBONE`, `_BASES`), which `tools/gen_chem_tables.py` checks.

Cut with `python3 -I tools/extract_arnott.py nucgen.dat`, which checks the SHA-256 of the copy read:

| Copy read | SHA-256 |
|---|---|
| https://github.com/ParmEd/ParmEd, `test/files/parm/nucgen.dat` (master at d3140213260b, 2026-11-26) | `2442dcb7769b00cec209c6a81710e9b26fad7e4a23e2b81389a23b56bad812bc` |
| https://github.com/choderalab/ambermini, `share/amber/dat/leap/parm/nucgen.dat` (master, 2026-11-26): identical | same |

The block lists two placeholder hydrogens on P's position (set `placeholder`, no PDB name), the
backbone (OA / OB = OP1 / OP2, O1' = O4', the old names) and then the bases A, G, C, U without
naming them; `tools/extract_arnott.py` assigns each run of rows to its base and checks it against the
atoms of the CCD's standard residue (`../ccd/ribonucleotide_atoms.tsv`). `tools/gen_chem_tables.py`
checks the OA / OB assignment: with OA as OP1 the chiral volume V(P; O5', OP1, OP2) is -2.46 A^3, the
sign of all 20 phosphates of PDB 4JRD (`../pdb`), which the wwPDB names by the PDB's OP1 / OP2 rule.
The table's base hydrogens are kept as facts but not used (the builder places hydrogens by rule).
