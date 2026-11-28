# AmberTools file and unit names of the Amber step types

`leap_files.tsv` records the names the Amber generators write into tleap / ParmEd inputs and build
scripts: the leaprc file of each RNA force field and water model, each water model's solvent box
unit, whether it has an extra point (four-site), its ParmEd `add12_6_4 watermodel` name and 12-6-4
ion frcmod, the residue and ion names of the strip masks, and the AmberTools files the embedded
helpers read at run time. Each row says where the name was found and the SHA-256 of that file; an
empty value records that AmberTools / ParmEd has none (checked: no such file or model exists).
Only names are recorded (facts); no AmberTools or ParmEd content is copied.

Cut with `python3 -I tools/extract_leap_files.py AMBERTOOLS PARMED`, which also checks per water
model that its leaprc names the water unit WAT and HOH, that the box unit is defined in a library
the leaprc loads, that its residues are named WAT and that it has atoms of type EP exactly when it
is a four-site model. `tools/gen_workflow_tables.py` builds `src/workflow/tables.inc` from it and
checks MDATK 0.0.4.1's tables (`mdatk/workflows/amber.py` FORCE_FIELDS, WATERS, FOUR_SITE_WATERS,
DEFAULT_STRIP, SOLVENT_MASK; "checked against AmberTools 26 dat/leap") equal it: they do.

| Source | Retrieved | SHA-256 | Licence |
|---|---|---|---|
| AmberTools 26.0, conda-forge build `linux-64/ambertools-26.0-cuda_None_nompi_py313h34d38f2_100.conda` (uploaded 2026-07-15), https://conda.anaconda.org/conda-forge/linux-64/ambertools-26.0-cuda_None_nompi_py313h34d38f2_100.conda; its `dat/leap` tree | 2026-11-28 | `218541561235f1425af60f1d5d75581ca4828d6609b53aa453a87550be95d507` (the package) | GPL-3.0-or-later AND LGPL-3.0-or-later AND BSD-3-Clause AND MIT (the package's metadata); only file and unit names are used (facts) |
| ParmEd 4.3.1, PyPI sdist `parmed-4.3.1.tar.gz`, https://files.pythonhosted.org/packages/a0/f5/65d3c930767aece106b97443eda3a37de8b9e6a183b7d999a670b50c161f/parmed-4.3.1.tar.gz; `parmed/tools/actions.py` (`add12_6_4`'s `_supported_wms`: TIP3P, TIP4PEW, SPCE, OPC3, OPC, FB3, FB4) | 2026-11-28 | `1f26bfcad2d9e1344fda343fe0199fd36ae9a02f64669010bf47e30f9ce0918c` (the sdist) | LGPL (PKG-INFO); only the water-model names are used (facts) |

Not AmberTools names, so not in the table: `TIP4P2005BOX`, `frcmod.tip4p2005`, `tip4p2005box.off`
(written by the tip4p2005 helper), `RNA_CP.lib`, `frcmod.cp`, `lj_1264_pol_cp.dat` (the cp-lib
helper), and the L-RNA files a user provides (`nucleic12_LRNA.lib`, `parm10_LRNA.dat` of Dudek &
Trylska; `frcmod.99bsc0-chiol3-CaseP`, the default of the own-CP-frcmod option), none of which
AmberTools 26 ships.
