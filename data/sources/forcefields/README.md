# Force-field residue and atom names

| File | Origin | Retrieved | Licence |
|---|---|---|---|
| `amber24_manual_names.tsv` | Amber 2024 Reference Manual, https://ambermd.org/doc12/Amber24.pdf: section 3.6 (solvents.lib water names), 12.2 (HID/HIE/HIP, CYX, ASH/GLH/LYN, ACE/NHE/NME), 13.4.2 (N/C-prefixed terminal amino acids), 13.4.3 (5/3/N nucleic end names) | 2026-11-24 | names only (facts); hand-typed from the manual text, no text copied |
| `gromacs_manual_names.tsv` | GROMACS reference manual, https://manual.gromacs.org/current/ : Table 12 "Internal GROMACS residue naming convention" (pdb2gmx input files) and the SOL molecule type of the topology examples | 2026-11-24 | names only (facts); no GROMACS data file (LGPL) was used |
| `openmm_residue_atoms.tsv` | OpenMM 8.1.2 force-field XML, https://github.com/openmm/openmm/tree/8.1.2/wrappers/python/openmm/app/data : `amber14/protein.ff14SB.xml`, `amber14/DNA.OL15.xml`, `amber14/RNA.OL3.xml`, `amber14/tip3p.xml`, `amber14/tip4pew.xml`, `amber14/opc.xml`, `tip5p.xml`, `charmm36.xml`, `charmm36/water.xml`, `charmm36/tip4pew.xml`, `charmm36/tip5p.xml`, `charmm36/spce.xml` | 2026-11-24 | MIT (OpenMM Python application layer, © Stanford University and the Authors) |
| `openmm_pdbnames_aliases.tsv` | OpenMM 8.1.2 `pdbNames.xml` (same directory): residue-name aliases and the nucleic and water atom-name aliases | 2026-11-24 | MIT |

`openmm_residue_atoms.tsv` keeps the residue names of the protein files, the heavy-atom names of the
nucleic residues, and every atom of the solvent/ion residues, each with the element of its atom
type ('' for virtual sites); cut with `python3 -I tools/extract_openmm.py DIR`. The amber14 files
are OpenMM's conversion of AmberTools' ff14SB, OL15, OL3 and water/ion libraries, the charmm36 files
its conversion of the CHARMM36 toppar (c36, Aug 2015); only names and elements are taken from them,
no parameters.

Not found in an admissible source and therefore not included: GROMACS's TIP4P virtual-site name
`MW` and a generic CHARMM lone-pair atom name `LP` (only `LP1`/`LP2` of TIP5P and `OM` of TIP4P-Ew
appear in the CHARMM36 files above), and the older Amber RNA names `RA`, `RU`, `RG`, `RC`.
