# Force-field residue and atom names

| File | Origin | Retrieved | Licence |
|---|---|---|---|
| `amber24_manual_names.tsv` | Amber 2024 Reference Manual, https://ambermd.org/doc12/Amber24.pdf: section 3.6 (solvents.lib water names), 12.2 (HID/HIE/HIP, CYX, ASH/GLH/LYN, ACE/NHE/NME), 13.4.2 (N/C-prefixed terminal amino acids), 13.4.3 (5/3/N nucleic end names) | 2026-11-24 | names only (facts); hand-typed from the manual text, no text copied |
| `gromacs_manual_names.tsv` | GROMACS reference manual, https://manual.gromacs.org/current/ : Table 12 "Internal GROMACS residue naming convention" (pdb2gmx input files) and the SOL molecule type of the topology examples | 2026-11-24 | names only (facts); no GROMACS data file (LGPL) was used |
| `openmm_residue_atoms.tsv` | OpenMM 8.1.2 force-field XML, https://github.com/openmm/openmm/tree/8.1.2/wrappers/python/openmm/app/data : `amber14/protein.ff14SB.xml`, `amber14/DNA.OL15.xml`, `amber14/RNA.OL3.xml`, `amber14/tip3p.xml`, `amber14/tip4pew.xml`, `amber14/opc.xml`, `tip5p.xml`, `charmm36.xml`, `charmm36/water.xml`, `charmm36/tip4pew.xml`, `charmm36/tip5p.xml`, `charmm36/spce.xml` | 2026-11-24 | MIT (OpenMM Python application layer, © Stanford University and the Authors) |
| `openmm_pdbnames_aliases.tsv` | OpenMM 8.1.2 `pdbNames.xml` (same directory): residue-name aliases and the nucleic and water atom-name aliases | 2026-11-24 | MIT |
| `openmm_amber14_templates.tsv` | OpenMM 8.1.2 `amber14/protein.ff14SB.xml`, `amber14/DNA.OL15.xml`, `amber14/RNA.OL3.xml`, `amber14/tip3p.xml` (same directory): every residue template, one row per atom: name, atom type, Amber class, element, partial charge, and its type's NonbondedForce sigma (nm) / epsilon (kJ/mol), as the files write them | 2026-11-25 | MIT |
| `openmm_pdbnames_residue_atom_aliases.tsv` | OpenMM 8.1.2 `pdbNames.xml`: the residue-specific atom-name aliases (ILE `CD` -> `CD1`, DT `C5M` -> `C7`, ...) | 2026-11-25 | MIT |
| `amber_lj_types.tsv` | numbers read from three AmberTools parameter files (below): R* (Angstrom) and epsilon (kcal/mol) of the Amber atom types the energy analysis' typed layer uses, with the file and line | 2026-11-25 | facts (numbers with their provenance; no file copied) |

`openmm_residue_atoms.tsv` keeps the residue names of the protein files, the heavy-atom names of the
nucleic residues, and every atom of the solvent/ion residues, each with the element of its atom
type ('' for virtual sites); cut with `python3 -I tools/extract_openmm.py DIR`. The amber14 files
are OpenMM's conversion of AmberTools' ff14SB, OL15, OL3 and water/ion libraries, the charmm36 files
its conversion of the CHARMM36 toppar (c36, Aug 2015); only names and elements are taken from them
for the name tables; the energy analysis' parameters (below) take charges and sigma / epsilon from
the four amber14 files.

Not found in an admissible source and therefore not included: GROMACS's TIP4P virtual-site name
`MW` and a generic CHARMM lone-pair atom name `LP` (only `LP1`/`LP2` of TIP5P and `OM` of TIP4P-Ew
appear in the CHARMM36 files above), and the older Amber RNA names `RA`, `RU`, `RG`, `RC`.

## Parameters of the energy analysis (docs/DIFFERENCES.md AN-25 - AN-27)

`openmm_amber14_templates.tsv` and `openmm_pdbnames_residue_atom_aliases.tsv` are cut by
`python3 -I tools/extract_openmm.py DIR` too. The four amber14 files are OpenMM's conversion of
AmberTools 17's `leaprc.protein.ff14SB`, `leaprc.DNA.OL15`, `leaprc.RNA.OL3` and the TIP3P water
with `frcmod.ionsjc_tip3p`, `frcmod.ions234lm_126_tip3p` and `atomic_ions.lib` (their `<Info>`
blocks). Files read (https://github.com/openmm/openmm/tree/8.1.2/wrappers/python/openmm/app/data,
fetched again on 2026-11-25 and identical to the 2026-11-24 copies), SHA-256:

| File | SHA-256 |
|---|---|
| `amber14/protein.ff14SB.xml` | `d9f9779c09d67cd5f8bc657692f174ffab14c469dfd06d560ac1899fa7e976b8` |
| `amber14/DNA.OL15.xml` | `200e7cd116481fb47ed50bfd29deefb84a1e488d73d262c3f3f278b36b80ffa5` |
| `amber14/RNA.OL3.xml` | `61b793a49ded3c9c3bd9fbab328bc1b20b2b4ce7bc2358d537892ca64b085127` |
| `amber14/tip3p.xml` | `3f4b188dbcb6c02863230eaca231e927fb6bf3307ce947d8a50d0f46f6dd83d9` |
| `pdbNames.xml` | `f36bc24f85d66c518a7cd3886699986ae4941c89023ceea562500d69ddfb16ee` |

`amber_lj_types.tsv` is cut by `python3 -I tools/extract_amber_lj.py DIR`, which checks these
SHA-256 sums; only the two numbers of each listed type are kept (the parameter files' licence is
not relied on: no file, comment or other content is copied, as MDATK's later generic tables were
made):

| File | Copy read | SHA-256 |
|---|---|---|
| `parm10.dat` (section `MOD4 RE`) | https://github.com/ParmEd/ParmEd, `test/files/parm/parm10.dat` (master, 2026-11-25) | `1ca292d5679e35450ea111a8d1b9442b5ece8518e824a8e43ed290dfff1a2ae6` |
| `frcmod.ionsjc_tip3p` (`NONBON`) | https://github.com/ParmEd/ParmEd, `test/files/parm/frcmod.ionsjc_tip3p` (master, 2026-11-25) | `eb7a1a59bf356581f89534287dd59d1bebbb2e5de7176d7c7dd69cb1d4928a55` |
| `frcmod.ions234lm_126_tip3p` (`NONBON`) | https://github.com/choderalab/ambermini, `share/amber/dat/leap/parm/` (master, 2026-11-25) | `f8bdb3df35ebdd9047b2e63786b4147265d5f73eb9b5197be3cb71027a2e25f3` |

These are the copies, and SHA-256 prefixes (`1ca292d5679e3545`, `eb7a1a59bf356581`,
`f8bdb3df35ebdd90`), that the documentation of MDATK's later generic Lennard-Jones tables (same
owner) records; every value it lists was checked against these lines. Types kept: `C`, `CT`, `N`,
`O`, `OH`, `OS`, `OW`, `P`, `S`, `F`, `Cl`, `Br`, `I`, `H`, `HO`, `HS`, `H1` (parm10.dat), `Li+`,
`Na+`, `K+`, `Rb+`, `Cs+`, `F-`, `Cl-`, `Br-`, `I-` (Joung & Cheatham 2008, TIP3P set), `Mg2+`,
`Ca2+`, `Sr2+`, `Ba2+`, `Mn2+`, `Fe2+`, `Co2+`, `Ni2+`, `Cu2+`, `Zn2+`, `Cd2+` (Li et al. 2013,
12-6 CM set for TIP3P). Cross-check: tip3p.xml's ions, converted back (R* = sigma 2^(1/6) / 2,
epsilon / 4.184), equal these rows to all printed digits.

`tools/gen_energy_lj.py` writes `src/analysis/energy_lj_tables.inc` from the three tables.
