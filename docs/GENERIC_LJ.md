# Lennard-Jones parameters of the energy analysis without parm7

The *Intermolecular energy* analysis needs Lennard-Jones parameters (and partial charges) for every atom. With
an Amber **parm7** topology and `analyses.energy.use_parm7` it uses that topology's own table, exactly
(`lj_source: parm7` in the result's metadata). Any other topology (PDB, GRO, PSF, a parm7 file that cannot be
read, or `use_parm7 = false`) gets its parameters in two layers:

1. **Residue templates.** A residue that matches a residue template of Amber's force fields takes every
   atom's parameters from that template: the force field's own values per atom type.
2. **Amber atom types.** Every other atom takes the values of one Amber atom type, chosen by its element and,
   where they are known, its bonded neighbours. This is the design of MDATK's generic table (MDATK 0.1.1,
   `docs/GENERIC_LJ.md` there, same owner), with the refinements listed below.

This file is the record of both layers: every number, the file it was read from (with its SHA-256), and the
publication behind it. The code is `src/analysis/energy_parameters.cpp` on the tables
`tools/gen_energy_lj.py` generates (`src/analysis/energy_lj_tables.inc`) from the extracts in
`data/sources/forcefields/` (`openmm_amber14_templates.tsv`, `amber_lj_types.tsv`; how they were cut:
`data/sources/forcefields/README.md`). The behaviour is described in `docs/DIFFERENCES.md` (AN-25 - AN-27,
CIT-04). The test `EnergyParameters.GenericLjDocument` (tests/analysis/test_energy.cpp) keeps the tables of this
file identical to those extracts.

The result's metadata says what was used: `lj_source` (`templates`, `typed` or `templates+typed`), `lj_atoms`
(atoms per layer), `lj_references` (the citation keys of the numbers used), `lj_default_elements` (elements
that fell back to carbon's values) and `lj_table`, MDATK 0.1.1's key, naming this table:
`Amber force-field residue templates (OpenMM amber14: ff14SB, OL15, OL3, TIP3P with Joung-Cheatham and Li-Merz
ions), then Amber parm10 atom types, Joung-Cheatham and Li-Merz (CM) ions for TIP3P (docs/GENERIC_LJ.md)`.
With parm7 parameters there is no `lj_table`.

## Units and combination

Every row is **R\*** = Rmin/2 in Å and the well depth **ε** in kcal/mol, as Amber's parameter files give
them. OpenMM's force-field files give σ in nm and ε in kJ/mol; they are converted as R\* = σ 2^(1/6) / 2 and
ε / 4.184, which gives Amber's values to about 1e-15 relative. A pair *i*, *j* at distance *r* has

    Rmin = R*_i + R*_j        ε_ij = sqrt(ε_i ε_j)        (Lorentz-Berthelot, as in Amber)
    E = ε_ij [ (Rmin / r)^12 - 2 (Rmin / r)^6 ]  =  A / r^12 - B / r^6,   A = ε_ij Rmin^12,   B = 2 ε_ij Rmin^6

Distances below 0.5 Å are taken as 0.5 Å. The electrostatic term is 332.0522 q_i q_j / (ε_r r) kcal/mol with
the topology's charges, or, for a topology without charges, the charges of the templates a residue matches
completely (AN-27).

## Layer 1: residue templates

The templates are those of OpenMM 8.1.2's conversion of AmberTools' force fields (`amber14/`):
`protein.ff14SB.xml` (ff14SB), `DNA.OL15.xml` (OL15), `RNA.OL3.xml` (OL3) and `tip3p.xml` (TIP3P water with the
Joung-Cheatham monovalent ions and the Li-Merz ions of `frcmod.ions234lm_126_tip3p` and the 3+ / 4+ set).

How a residue gets a template:

1. **Candidates:** the template named like the residue (upper case), with its terminal variants when it is the
   middle form (Amber's N- / C-prefixed amino acids, the 5 / 3 / N-suffixed nucleotides); a residue name that
   is no template stands for the templates `pdbNames.xml` makes it an alias of (HIS, HSD -> HID / HIE / HIP;
   WAT, SOL, TIP3 -> HOH; ADE -> A ...); a one-atom residue is an ion template named like the residue or the
   atom (as written, upper case, without charge signs and digits: `Na+` -> NA) of the same element.
2. **Atoms:** an atom maps to the template atom of its name, else by `pdbNames.xml`'s aliases (residue-specific
   first, then the Protein / Nucleic / HOH ones: `O1P` -> OP1, `HO5'` <- `H5T` ...), one template atom per atom;
   its element, when known, must be the template atom's. A residue *matches* a template when every atom maps,
   *completely* when every template atom is mapped too (hydrogens included). Extra points are left out.
3. **Choice:** exactly one complete match: that template, whatever the residue's position. Otherwise the
   templates that contain the residue are narrowed by its position in its chain when it has terminal variants
   (start -> N / 5', end -> C / 3', both -> the nucleic N form, neither -> the middle form). A residue starts a
   chain when it does not continue the residue before it (by index), ends one when the next does not continue
   it. Residue r + 1 continues r unless the evidence agrees on a break: a new chain (another segment, or
   another chain ID where the topology has chain IDs), or a gap or reset in the numbering *and* no bond between
   the two (the topology's bonds where either residue has any, else bonds guessed on the first analysed frame).
4. If the remaining templates give every atom the same R\* and ε (HID / HIE / HIP without hydrogens do), the
   residue takes them; else (none left, or they disagree) its atoms go to layer 2.

The templates use these Lennard-Jones classes. Every protein and nucleic class has the R\* / ε of a type of
`parm94.dat` (checked class by class against ParmEd's copy, SHA-256 `34021e968833329f...`; parm99 and parm10
kept the values), so these numbers are cited to Cornell et al. 1995 (CIT-04). A class whose ε is 0 has no
Lennard-Jones term (OpenMM writes σ = 1 Å for it; "-" below).

### Table 1: classes of the protein and nucleic templates

| Class | R\* (Å) | ε (kcal/mol) | Files | Reference |
|---|---|---|---|---|
| `2C` | 1.908 | 0.1094 | `amber14/protein.ff14SB.xml` | Cornell 1995 |
| `3C` | 1.908 | 0.1094 | `amber14/protein.ff14SB.xml` | Cornell 1995 |
| `C` | 1.908 | 0.086 | `amber14/DNA.OL15.xml`, `amber14/RNA.OL3.xml`, `amber14/protein.ff14SB.xml` | Cornell 1995 |
| `C*` | 1.908 | 0.086 | `amber14/protein.ff14SB.xml` | Cornell 1995 |
| `C1` | 1.908 | 0.086 | `amber14/DNA.OL15.xml` | Cornell 1995 |
| `C2` | 1.908 | 0.086 | `amber14/DNA.OL15.xml` | Cornell 1995 |
| `C4` | 1.908 | 0.086 | `amber14/RNA.OL3.xml` | Cornell 1995 |
| `C5` | 1.908 | 0.086 | `amber14/RNA.OL3.xml` | Cornell 1995 |
| `C7` | 1.908 | 0.1094 | `amber14/DNA.OL15.xml` | Cornell 1995 |
| `C8` | 1.908 | 0.1094 | `amber14/protein.ff14SB.xml` | Cornell 1995 |
| `CA` | 1.908 | 0.086 | `amber14/DNA.OL15.xml`, `amber14/RNA.OL3.xml`, `amber14/protein.ff14SB.xml` | Cornell 1995 |
| `CB` | 1.908 | 0.086 | `amber14/DNA.OL15.xml`, `amber14/RNA.OL3.xml`, `amber14/protein.ff14SB.xml` | Cornell 1995 |
| `CC` | 1.908 | 0.086 | `amber14/protein.ff14SB.xml` | Cornell 1995 |
| `CI` | 1.908 | 0.1094 | `amber14/RNA.OL3.xml` | Cornell 1995 |
| `CJ` | 1.908 | 0.1094 | `amber14/DNA.OL15.xml` | Cornell 1995 |
| `CK` | 1.908 | 0.086 | `amber14/DNA.OL15.xml` | Cornell 1995 |
| `CM` | 1.908 | 0.086 | `amber14/DNA.OL15.xml` | Cornell 1995 |
| `CN` | 1.908 | 0.086 | `amber14/protein.ff14SB.xml` | Cornell 1995 |
| `CO` | 1.908 | 0.086 | `amber14/protein.ff14SB.xml` | Cornell 1995 |
| `CP` | 1.908 | 0.086 | `amber14/RNA.OL3.xml` | Cornell 1995 |
| `CQ` | 1.908 | 0.086 | `amber14/DNA.OL15.xml`, `amber14/RNA.OL3.xml` | Cornell 1995 |
| `CR` | 1.908 | 0.086 | `amber14/protein.ff14SB.xml` | Cornell 1995 |
| `CS` | 1.908 | 0.086 | `amber14/RNA.OL3.xml` | Cornell 1995 |
| `CT` | 1.908 | 0.1094 | `amber14/DNA.OL15.xml`, `amber14/RNA.OL3.xml`, `amber14/protein.ff14SB.xml` | Cornell 1995 |
| `CV` | 1.908 | 0.086 | `amber14/protein.ff14SB.xml` | Cornell 1995 |
| `CW` | 1.908 | 0.086 | `amber14/protein.ff14SB.xml` | Cornell 1995 |
| `CX` | 1.908 | 0.1094 | `amber14/protein.ff14SB.xml` | Cornell 1995 |
| `H` | 0.6 | 0.0157 | `amber14/DNA.OL15.xml`, `amber14/RNA.OL3.xml`, `amber14/protein.ff14SB.xml` | Cornell 1995 |
| `H1` | 1.387 | 0.0157 | `amber14/DNA.OL15.xml`, `amber14/RNA.OL3.xml`, `amber14/protein.ff14SB.xml` | Cornell 1995 |
| `H2` | 1.287 | 0.0157 | `amber14/DNA.OL15.xml`, `amber14/RNA.OL3.xml` | Cornell 1995 |
| `H4` | 1.409 | 0.015 | `amber14/DNA.OL15.xml`, `amber14/RNA.OL3.xml`, `amber14/protein.ff14SB.xml` | Cornell 1995 |
| `H5` | 1.359 | 0.015 | `amber14/DNA.OL15.xml`, `amber14/RNA.OL3.xml`, `amber14/protein.ff14SB.xml` | Cornell 1995 |
| `HA` | 1.459 | 0.015 | `amber14/DNA.OL15.xml`, `amber14/RNA.OL3.xml`, `amber14/protein.ff14SB.xml` | Cornell 1995 |
| `HC` | 1.487 | 0.0157 | `amber14/DNA.OL15.xml`, `amber14/protein.ff14SB.xml` | Cornell 1995 |
| `HO` | - | 0 | `amber14/DNA.OL15.xml`, `amber14/RNA.OL3.xml`, `amber14/protein.ff14SB.xml` | Cornell 1995 |
| `HP` | 1.1 | 0.0157 | `amber14/protein.ff14SB.xml` | Cornell 1995 |
| `HS` | 0.6 | 0.0157 | `amber14/protein.ff14SB.xml` | Cornell 1995 |
| `N` | 1.824 | 0.17 | `amber14/protein.ff14SB.xml` | Cornell 1995 |
| `N*` | 1.824 | 0.17 | `amber14/DNA.OL15.xml`, `amber14/RNA.OL3.xml` | Cornell 1995 |
| `N2` | 1.824 | 0.17 | `amber14/DNA.OL15.xml`, `amber14/RNA.OL3.xml`, `amber14/protein.ff14SB.xml` | Cornell 1995 |
| `N3` | 1.824 | 0.17 | `amber14/protein.ff14SB.xml` | Cornell 1995 |
| `NA` | 1.824 | 0.17 | `amber14/DNA.OL15.xml`, `amber14/RNA.OL3.xml`, `amber14/protein.ff14SB.xml` | Cornell 1995 |
| `NB` | 1.824 | 0.17 | `amber14/DNA.OL15.xml`, `amber14/RNA.OL3.xml`, `amber14/protein.ff14SB.xml` | Cornell 1995 |
| `NC` | 1.824 | 0.17 | `amber14/DNA.OL15.xml`, `amber14/RNA.OL3.xml` | Cornell 1995 |
| `O` | 1.6612 | 0.21 | `amber14/DNA.OL15.xml`, `amber14/RNA.OL3.xml`, `amber14/protein.ff14SB.xml` | Cornell 1995 |
| `O2` | 1.6612 | 0.21 | `amber14/DNA.OL15.xml`, `amber14/RNA.OL3.xml`, `amber14/protein.ff14SB.xml` | Cornell 1995 |
| `OH` | 1.721 | 0.2104 | `amber14/DNA.OL15.xml`, `amber14/RNA.OL3.xml`, `amber14/protein.ff14SB.xml` | Cornell 1995 |
| `OS` | 1.6837 | 0.17 | `amber14/DNA.OL15.xml`, `amber14/RNA.OL3.xml` | Cornell 1995 |
| `P` | 2.1 | 0.2 | `amber14/DNA.OL15.xml`, `amber14/RNA.OL3.xml` | Cornell 1995 |
| `S` | 2 | 0.25 | `amber14/protein.ff14SB.xml` | Cornell 1995 |
| `SH` | 2 | 0.25 | `amber14/protein.ff14SB.xml` | Cornell 1995 |

### Table 2: TIP3P water

| Atom | R\* (Å) | ε (kcal/mol) | File | Reference |
|---|---|---|---|---|
| `H` | - | 0 | `amber14/tip3p.xml` | Jorgensen 1983 |
| `O` | 1.7683 | 0.152 | `amber14/tip3p.xml` | Jorgensen 1983 |

### Table 3: ion templates of `tip3p.xml`

Monovalent ions: Joung & Cheatham 2008 (TIP3P set, `frcmod.ionsjc_tip3p`); divalent: Li et al. 2013 (12-6
"CM" set for TIP3P, `frcmod.ions234lm_126_tip3p`); 3+ and 4+: Li et al. 2015 (12-6 set for TIP3P).

| Ion | R\* (Å) | ε (kcal/mol) | File | Reference |
|---|---|---|---|---|
| `Cs+` | 1.976 | 0.4065394 | `amber14/tip3p.xml` | Joung 2008 |
| `K+` | 1.705 | 0.1936829 | `amber14/tip3p.xml` | Joung 2008 |
| `Li+` | 1.025 | 0.0279896 | `amber14/tip3p.xml` | Joung 2008 |
| `Na+` | 1.369 | 0.0874393 | `amber14/tip3p.xml` | Joung 2008 |
| `Rb+` | 1.813 | 0.3278219 | `amber14/tip3p.xml` | Joung 2008 |
| `Br-` | 2.608 | 0.0586554 | `amber14/tip3p.xml` | Joung 2008 |
| `Cl-` | 2.513 | 0.035591 | `amber14/tip3p.xml` | Joung 2008 |
| `F-` | 2.303 | 0.003364 | `amber14/tip3p.xml` | Joung 2008 |
| `I-` | 2.86 | 0.0536816 | `amber14/tip3p.xml` | Joung 2008 |
| `Ag2+` | 1.336 | 0.00770969 | `amber14/tip3p.xml` | Li 2013 |
| `Ba2+` | 2.019 | 0.40664608 | `amber14/tip3p.xml` | Li 2013 |
| `Be2+` | 0.956 | 3.95e-06 | `amber14/tip3p.xml` | Li 2013 |
| `Ca2+` | 1.649 | 0.1059287 | `amber14/tip3p.xml` | Li 2013 |
| `Cd2+` | 1.412 | 0.01773416 | `amber14/tip3p.xml` | Li 2013 |
| `Co2+` | 1.299 | 0.00483892 | `amber14/tip3p.xml` | Li 2013 |
| `Cr2+` | 1.346 | 0.00868178 | `amber14/tip3p.xml` | Li 2013 |
| `Cu2+` | 1.218 | 0.00148497 | `amber14/tip3p.xml` | Li 2013 |
| `Eu2+` | 1.802 | 0.21475916 | `amber14/tip3p.xml` | Li 2013 |
| `Fe2+` | 1.353 | 0.00941798 | `amber14/tip3p.xml` | Li 2013 |
| `Hg2+` | 1.407 | 0.0168671 | `amber14/tip3p.xml` | Li 2013 |
| `Mg2+` | 1.36 | 0.01020237 | `amber14/tip3p.xml` | Li 2013 |
| `Mn2+` | 1.407 | 0.0168671 | `amber14/tip3p.xml` | Li 2013 |
| `Ni2+` | 1.255 | 0.0026232 | `amber14/tip3p.xml` | Li 2013 |
| `Pb2+` | 1.745 | 0.17018074 | `amber14/tip3p.xml` | Li 2013 |
| `Pd2+` | 1.303 | 0.00509941 | `amber14/tip3p.xml` | Li 2013 |
| `Pt2+` | 1.266 | 0.00307642 | `amber14/tip3p.xml` | Li 2013 |
| `Ra2+` | 2.019 | 0.40664608 | `amber14/tip3p.xml` | Li 2013 |
| `Sm2+` | 1.819 | 0.22878796 | `amber14/tip3p.xml` | Li 2013 |
| `Sn2+` | 1.666 | 0.11617738 | `amber14/tip3p.xml` | Li 2013 |
| `Sr2+` | 1.81 | 0.22132374 | `amber14/tip3p.xml` | Li 2013 |
| `V2+` | 1.364 | 0.01067299 | `amber14/tip3p.xml` | Li 2013 |
| `Yb2+` | 1.642 | 0.10185975 | `amber14/tip3p.xml` | Li 2013 |
| `Zn2+` | 1.271 | 0.00330286 | `amber14/tip3p.xml` | Li 2013 |
| `Al3+` | 1.369 | 0.01128487 | `amber14/tip3p.xml` | Li 2015 |
| `Ce3+` | 1.782 | 0.19865859 | `amber14/tip3p.xml` | Li 2015 |
| `Cr3+` | 1.415 | 0.01827024 | `amber14/tip3p.xml` | Li 2015 |
| `Dy3+` | 1.637 | 0.09900804 | `amber14/tip3p.xml` | Li 2015 |
| `Er3+` | 1.635 | 0.09788018 | `amber14/tip3p.xml` | Li 2015 |
| `Eu3+` | 1.716 | 0.14920231 | `amber14/tip3p.xml` | Li 2015 |
| `Fe3+` | 1.443 | 0.02387506 | `amber14/tip3p.xml` | Li 2015 |
| `Gd3+` | 1.658 | 0.11129023 | `amber14/tip3p.xml` | Li 2015 |
| `In3+` | 1.491 | 0.03625449 | `amber14/tip3p.xml` | Li 2015 |
| `La3+` | 1.758 | 0.1799796 | `amber14/tip3p.xml` | Li 2015 |
| `Lu3+` | 1.625 | 0.09235154 | `amber14/tip3p.xml` | Li 2015 |
| `Nd3+` | 1.724 | 0.15486311 | `amber14/tip3p.xml` | Li 2015 |
| `Pr3+` | 1.78 | 0.19707431 | `amber14/tip3p.xml` | Li 2015 |
| `Sm3+` | 1.711 | 0.14571499 | `amber14/tip3p.xml` | Li 2015 |
| `Tb3+` | 1.671 | 0.11928915 | `amber14/tip3p.xml` | Li 2015 |
| `Tl3+` | 1.571 | 0.0657303 | `amber14/tip3p.xml` | Li 2015 |
| `Tm3+` | 1.647 | 0.10475707 | `amber14/tip3p.xml` | Li 2015 |
| `Y3+` | 1.63 | 0.09509276 | `amber14/tip3p.xml` | Li 2015 |
| `Ce4+` | 1.766 | 0.18612361 | `amber14/tip3p.xml` | Li 2015 |
| `Hf4+` | 1.6 | 0.07934493 | `amber14/tip3p.xml` | Li 2015 |
| `Pu4+` | 1.752 | 0.17542802 | `amber14/tip3p.xml` | Li 2015 |
| `Th4+` | 1.77 | 0.18922704 | `amber14/tip3p.xml` | Li 2015 |
| `U4+` | 1.792 | 0.20665151 | `amber14/tip3p.xml` | Li 2015 |
| `Zr4+` | 1.609 | 0.0838924 | `amber14/tip3p.xml` | Li 2015 |

## Layer 2: Amber atom types

An atom of a residue no template matches (ligands, modified residues, ions with other names, a 5'-terminal
nucleotide that keeps its phosphate) takes the row of one Amber atom type:

1. **A hydrogen** takes the type of the atom it is bonded to: on N `H`, on O `HO`, on S `HS`, otherwise `H1`
   (table 4). The bonded atom comes from the topology's bonds; without them it is the nearest non-hydrogen atom
   of the hydrogen's residue within 1.45 Å on the first analysed frame; with none, the hydrogen is taken as on
   nitrogen (the small radius, which cannot turn a contact repulsive).
2. **An atom alone in its residue** whose element has a TIP3P ion set (table 5) is that ion.
3. **C, N, O, P, S, F, Cl, Br, I** take their parm10 type of table 4, refined by the bonded neighbours where
   they are known (the topology's bonds where the residue has any, else bonds guessed on the first analysed
   frame): a carbon with four neighbours `CT` (sp3, as Amber types it) instead of `C`; an oxygen with two
   neighbours `OS` (ether, ester, phosphate ester), `OH` when one of them is a hydrogen, `OW` (TIP3P's water
   oxygen) when both are. Without hydrogens a carbon rarely has four neighbours, so it stays `C`. MDATK 0.1.1
   gives every carbon `C` and every oxygen `O`.
4. **A metal of table 5 inside a larger residue** takes its ion's row.
5. **Extra points** (virtual sites, lone pairs: an atom whose element is no chemical element and whose name is
   a water virtual site of OpenMM's force fields (M, M1, M2, LP1, LP2, OM), GROMACS' TIP4P `MW`, or starts with
   `EP` or `LP`) take `EP`: no Lennard-Jones term. MDATK 0.1.1 gives them carbon's row.
6. **Any other element** takes carbon's row (`C`), with the warning `energy: no Lennard-Jones parameters for
   element(s) 'X' (N atoms); carbon's (Amber type C) are used` and `lj_default_elements` in the metadata.

The element is the topology's, or is guessed from the atom name where the topology has none.

### Table 4: atoms of molecules (`parm10.dat`, section `MOD4 RE`)

`parm10.dat` is the parameter file of ff10, ff14SB, OL3 and OL15. The numbers are written as the file writes
them; the last column but one is the line of the file they were read from.

| Type | R\* (Å) | ε (kcal/mol) | File | Line | Reference |
|---|---|---|---|---|---|
| `C` | 1.9080 | 0.0860 | `parm10.dat` | 987 | Cornell 1995 |
| `CT` | 1.9080 | 0.1094 | `parm10.dat` | 985 | Cornell 1995 |
| `N` | 1.8240 | 0.1700 | `parm10.dat` | 988 | Cornell 1995 |
| `O` | 1.6612 | 0.2100 | `parm10.dat` | 975 | Cornell 1995 |
| `OH` | 1.7210 | 0.2104 | `parm10.dat` | 978 | Cornell 1995 |
| `OS` | 1.6837 | 0.1700 | `parm10.dat` | 979 | Cornell 1995 |
| `OW` | 1.7683 | 0.1520 | `parm10.dat` | 977 | Jorgensen 1983 |
| `P` | 2.1000 | 0.2000 | `parm10.dat` | 992 | Cornell 1995 |
| `S` | 2.0000 | 0.2500 | `parm10.dat` | 990 | Cornell 1995 |
| `F` | 1.75 | 0.061 | `parm10.dat` | 996 | Cornell 1995 |
| `Cl` | 1.948 | 0.265 | `parm10.dat` | 997 | Wang 2000 |
| `Br` | 2.22 | 0.320 | `parm10.dat` | 998 | Wang 2000 |
| `I` | 2.35 | 0.40 | `parm10.dat` | 999 | Cornell 1995 |
| `H` | 0.6000 | 0.0157 | `parm10.dat` | 962 | Cornell 1995 |
| `HO` | 0.0000 | 0.0000 | `parm10.dat` | 963 | Cornell 1995 |
| `HS` | 0.6000 | 0.0157 | `parm10.dat` | 964 | Cornell 1995 |
| `H1` | 1.3870 | 0.0157 | `parm10.dat` | 966 | Cornell 1995 |
| `EP` | 0.00 | 0.0000 | `parm10.dat` | 1000 | - |

### Table 5: monatomic ions

What AmberTools' `leaprc.water.tip3p` loads: `frcmod.ionsjc_tip3p` for the monovalent ions and
`frcmod.ions234lm_126_tip3p` for the divalent ones (the 12-6 CM set). The water model is not known without a
parm7 topology, and ion parameters are fitted per water model: the TIP3P sets are used throughout (the same
publication's Na+ is 1.212 / 0.3526418 for SPC/E). Iron is taken as Fe2+ and copper as Cu2+.

| Ion | R\* (Å) | ε (kcal/mol) | File | Line | Reference |
|---|---|---|---|---|---|
| `Li+` | 1.025 | 0.0279896 | `frcmod.ionsjc_tip3p` | 14 | Joung 2008 |
| `Na+` | 1.369 | 0.0874393 | `frcmod.ionsjc_tip3p` | 15 | Joung 2008 |
| `K+` | 1.705 | 0.1936829 | `frcmod.ionsjc_tip3p` | 16 | Joung 2008 |
| `Rb+` | 1.813 | 0.3278219 | `frcmod.ionsjc_tip3p` | 17 | Joung 2008 |
| `Cs+` | 1.976 | 0.4065394 | `frcmod.ionsjc_tip3p` | 18 | Joung 2008 |
| `F-` | 2.303 | 0.0033640 | `frcmod.ionsjc_tip3p` | 19 | Joung 2008 |
| `Cl-` | 2.513 | 0.0355910 | `frcmod.ionsjc_tip3p` | 20 | Joung 2008 |
| `Br-` | 2.608 | 0.0586554 | `frcmod.ionsjc_tip3p` | 21 | Joung 2008 |
| `I-` | 2.860 | 0.0536816 | `frcmod.ionsjc_tip3p` | 22 | Joung 2008 |
| `Mg2+` | 1.360 | 0.01020237 | `frcmod.ions234lm_126_tip3p` | 63 | Li 2013 |
| `Ca2+` | 1.649 | 0.10592870 | `frcmod.ions234lm_126_tip3p` | 69 | Li 2013 |
| `Sr2+` | 1.810 | 0.22132374 | `frcmod.ions234lm_126_tip3p` | 73 | Li 2013 |
| `Ba2+` | 2.019 | 0.40664608 | `frcmod.ions234lm_126_tip3p` | 75 | Li 2013 |
| `Mn2+` | 1.407 | 0.01686710 | `frcmod.ions234lm_126_tip3p` | 65 | Li 2013 |
| `Fe2+` | 1.353 | 0.00941798 | `frcmod.ions234lm_126_tip3p` | 62 | Li 2013 |
| `Co2+` | 1.299 | 0.00483892 | `frcmod.ions234lm_126_tip3p` | 58 | Li 2013 |
| `Ni2+` | 1.255 | 0.00262320 | `frcmod.ions234lm_126_tip3p` | 55 | Li 2013 |
| `Cu2+` | 1.218 | 0.00148497 | `frcmod.ions234lm_126_tip3p` | 54 | Li 2013 |
| `Zn2+` | 1.271 | 0.00330286 | `frcmod.ions234lm_126_tip3p` | 57 | Li 2013 |
| `Cd2+` | 1.412 | 0.01773416 | `frcmod.ions234lm_126_tip3p` | 67 | Li 2013 |

## Where the numbers were read

OpenMM 8.1.2's force-field files, https://github.com/openmm/openmm/tree/8.1.2/wrappers/python/openmm/app/data
(MIT; read on 2026-11-24 and again on 2026-11-25, identical):

| File | SHA-256 |
|---|---|
| `amber14/protein.ff14SB.xml` | `d9f9779c09d67cd5f8bc657692f174ffab14c469dfd06d560ac1899fa7e976b8` |
| `amber14/DNA.OL15.xml` | `200e7cd116481fb47ed50bfd29deefb84a1e488d73d262c3f3f278b36b80ffa5` |
| `amber14/RNA.OL3.xml` | `61b793a49ded3c9c3bd9fbab328bc1b20b2b4ce7bc2358d537892ca64b085127` |
| `amber14/tip3p.xml` | `3f4b188dbcb6c02863230eaca231e927fb6bf3307ce947d8a50d0f46f6dd83d9` |
| `pdbNames.xml` | `f36bc24f85d66c518a7cd3886699986ae4941c89023ceea562500d69ddfb16ee` |

The three AmberTools parameter files of layer 2, read on 2026-11-25 from public copies (only the two numbers
of each listed type are kept; `tools/extract_amber_lj.py` checks the sums). They are the copies, and SHA-256
prefixes, that MDATK 0.1.1's `docs/GENERIC_LJ.md` records, and every value it lists is the same here:

| File | Copy read | SHA-256 |
|---|---|---|
| `parm10.dat` | github.com/ParmEd/ParmEd, `test/files/parm/parm10.dat` | `1ca292d5679e35450ea111a8d1b9442b5ece8518e824a8e43ed290dfff1a2ae6` |
| `frcmod.ionsjc_tip3p` | github.com/ParmEd/ParmEd, `test/files/parm/frcmod.ionsjc_tip3p` | `eb7a1a59bf356581f89534287dd59d1bebbb2e5de7176d7c7dd69cb1d4928a55` |
| `frcmod.ions234lm_126_tip3p` | github.com/choderalab/ambermini, `share/amber/dat/leap/parm/` | `f8bdb3df35ebdd9047b2e63786b4147265d5f73eb9b5197be3cb71027a2e25f3` |

Cross-check: the ions of `tip3p.xml` (table 3), converted back from σ and ε, equal the rows of table 5 to all
printed digits.

## Measured

MDATK's demo system (a 12-base-pair RNA duplex with hydrogens, strand against strand, first frame, 12 Å
cutoff), van der Waals energy in kcal/mol: -34.67 with the parm7 file's own parameters; -34.62 from the
templates (layer 1); -34.74 from the atom types alone (layer 2, every residue renamed); MDATK 0.1.1's table
-34.1; MDATK up to 0.1.0 (one row per element, every hydrogen `HP`) +145.4.

## References

The analysis cites, in its log line and `lj_references`, only the publications behind the numbers a run used
(docs/DIFFERENCES.md CIT-04). MDATK 0.1.1's citation resource `generic_lj` is known to MDBridge (keys
`cornell1995`, `wang2000`, `jorgensen1983`, `joung2008`, `li2013`: the references of layer 2) so that its
params logs resolve, but is never a project resource.

- **Cornell 1995** (the protein and nucleic templates' charges and Lennard-Jones classes; the parm94 types of
  table 4): Cornell W. D., Cieplak P., Bayly C. I., Gould I. R., Merz K. M., Ferguson D. M., Spellmeyer D. C.,
  Fox T., Caldwell J. W., Kollman P. A. *J. Am. Chem. Soc.* 1995, 117, 5179-5197. doi:10.1021/ja00124a002
- **Wang 2000** (`Cl`, `Br`, parm99): Wang J., Cieplak P., Kollman P. A. *J. Comput. Chem.* 2000, 21,
  1049-1074. doi:10.1002/1096-987X(200009)21:12<1049::AID-JCC3>3.0.CO;2-F
- **Jorgensen 1983** (TIP3P water, `OW`): Jorgensen W. L., Chandrasekhar J., Madura J. D., Impey R. W., Klein
  M. L. *J. Chem. Phys.* 1983, 79, 926-935. doi:10.1063/1.445869
- **Joung 2008** (monovalent ions): Joung I. S., Cheatham T. E. *J. Phys. Chem. B* 2008, 112, 9020-9041.
  doi:10.1021/jp8001614
- **Li 2013** (divalent ions): Li P., Roberts B. P., Chakravorty D. K., Merz K. M. *J. Chem. Theory Comput.*
  2013, 9, 2733-2748. doi:10.1021/ct400146w
- **Li 2015** (3+ and 4+ ions): Li P., Song L. F., Merz K. M. *J. Phys. Chem. B* 2015, 119, 883-895.
  doi:10.1021/jp505875v

Not traced further: the charges of the residues Amber added after ff94 (ASH, GLH, LYN, CYM, HYP), cited as
ff94's; the per-line comments of `parm10.dat` (OPLS, Veenstra et al., ...), which name the origin of single
rows as the Amber developers recorded it.
