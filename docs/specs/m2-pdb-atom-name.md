# PDB writer: placing the atom name in columns 13–16 — behaviour specification

Reference behaviour: MDAnalysis 2.11.0 PDB writer. Every example in the
table was produced by writing atoms with the given (residue name, atom name)
using MDAnalysis and reading back columns 13–16 of the ATOM record.

## 1. Inputs and output

* Input: the atom name and the residue name, exactly as stored. They are not
  stripped and not upper-cased, and comparisons are case-sensitive. The full
  residue name is used for the decision, even though the writer itself
  prints only its first 4 characters in columns 18–21.
* Output: a 4-character field for columns 13–16. It is either
  * **"column 13" alignment**: the name left-justified in the 4 columns and
    padded with spaces on the right, or
  * **"column 14" alignment**: one space, then the name left-justified in the
    remaining 3 columns and padded on the right.

## 2. Rules

Let *n* be the atom name, *r* the residue name, and |n| the number of
characters in n.

1. **Empty name**: four spaces.
2. **|n| ≥ 4**: the first 4 characters of n, unchanged (column 13). Nothing
   else is checked, so trailing spaces in n are kept (`C14 ` → `C14 `).
3. **|n| = 1**: column 14, giving space + n + two spaces (e.g. ` N  `). This
   holds even if r equals n, so residue `K` with atom `K` gives ` K  `.
4. **|n| = 2 or 3**: **column 13** if at least one condition in group A
   holds **and** the (r, n) pair is not in the exclude list (§3.5).
   Otherwise **column 14**.
   Group A:
   * A1. r equals n exactly;
   * A2. the first two characters of n are in the ion list (§3.1);
   * A3. n is exactly `UNK`;
   * A4. r is in the mercury-residue list (§3.2) and n starts with `HG`;
   * A5. r is in the chlorine-residue list (§3.3) and n starts with `CL`;
   * A6. the pair (r, n) is in the include list (§3.4).

## 3. Data lists

### 3.1 Ion prefixes (first two characters of the atom name)

```
FE AS ZN MG MN CO BR CU TA MO AL BE SE PT EU NI IR RH AU GD RU
```

`CA`, `NA`, `CL`, `K`, `HG`, `PD` and `CD` are **not** on this list.

### 3.2 Mercury residues (with an atom name starting `HG`)

```
CMH EMC MBO MMC HGB BE7 PMB
```

### 3.3 Chlorine residues (with an atom name starting `CL`)

```
0QE CPT DCE EAA IMN OCZ OMY OMZ UN9 1N1 2T8 393 3MY BMU CLM CP6 DB8 DIF
EFZ LUR RDC UCL XMM HLT IRE LCP PCI VGH
```

### 3.4 Include pairs (residue, atom) → column 13

```
(OEC, CA1)  (PLL, PD)  (OEX, CA1)
```

### 3.5 Exclude pairs (residue, atom) → always column 14 (checked last)

```
(C14, C14) (C15, C15) (F9F, F9F) (OAN, OAN) (BLM, NI) (BZG, CO)
(BZG, NI)  (VNL, CO1) (VNL, CO2) (PF5, FE1) (PF5, FE2) (UNL, UNL)
```

## 4. Examples

In the field column, `·` shows a space. "Col" is the column where the first
non-space character lands; "blank" means the field is all spaces.

| residue | atom | columns 13–16 | col | rule |
|---|---|---|---|---|
| ALA | N | `·N··` | 14 | 1 char |
| ALA | CA | `·CA·` | 14 | no A-condition |
| ALA | C | `·C··` | 14 | |
| ALA | O | `·O··` | 14 | |
| ALA | CB | `·CB·` | 14 | |
| ALA | HB1 | `·HB1` | 14 | |
| ALA | 1HB1 | `1HB1` | 13 | 4 chars |
| ALA | HB11 | `HB11` | 13 | 4 chars |
| LYS | HZ1 | `·HZ1` | 14 | |
| ALA | OXT | `·OXT` | 14 | |
| CA | CA | `CA··` | 13 | A1 |
| NA | NA | `NA··` | 13 | A1 |
| ALA | NA | `·NA·` | 14 | NA not an ion prefix |
| SOD | SOD | `SOD·` | 13 | A1 |
| CL | CL | `CL··` | 13 | A1 |
| ALA | CL | `·CL·` | 14 | |
| CLA | CLA | `CLA·` | 13 | A1 |
| ZN | ZN | `ZN··` | 13 | A1/A2 |
| ALA | ZN | `ZN··` | 13 | A2 |
| ALA | ZN1 | `ZN1·` | 13 | A2 |
| ALA | FE | `FE··` | 13 | A2 |
| HEM | FE | `FE··` | 13 | A2 |
| HEM | FE1 | `FE1·` | 13 | A2 |
| HEM | NA | `·NA·` | 14 | |
| HEM | CHA | `·CHA` | 14 | |
| MSE | SE | `SE··` | 13 | A2 |
| ALA | SE | `SE··` | 13 | A2 |
| ALA | MG | `MG··` | 13 | A2 |
| ALA | MN | `MN··` | 13 | A2 |
| ALA | CO | `CO··` | 13 | A2 |
| ALA | CO1 | `CO1·` | 13 | A2 |
| ALA | CU | `CU··` | 13 | A2 |
| ALA | BR | `BR··` | 13 | A2 |
| ALA | BRX | `BRX·` | 13 | A2 |
| ALA | AS | `AS··` | 13 | A2 |
| ALA | TA | `TA··` | 13 | A2 |
| ALA | MO | `MO··` | 13 | A2 |
| ALA | AL | `AL··` | 13 | A2 |
| ALA | BE | `BE··` | 13 | A2 |
| ALA | PT | `PT··` | 13 | A2 |
| ALA | EU | `EU··` | 13 | A2 |
| ALA | NI | `NI··` | 13 | A2 |
| ALA | IR | `IR··` | 13 | A2 |
| ALA | RH | `RH··` | 13 | A2 |
| ALA | AU | `AU··` | 13 | A2 |
| ALA | AUX | `AUX·` | 13 | A2 |
| ALA | GD | `GD··` | 13 | A2 |
| ALA | RU | `RU··` | 13 | A2 |
| ALA | FEX | `FEX·` | 13 | A2 |
| ALA | MGX | `MGX·` | 13 | A2 |
| ALA | fe | `·fe·` | 14 | case-sensitive |
| ALA | Fe | `·Fe·` | 14 | case-sensitive |
| fe | fe | `fe··` | 13 | A1 |
| ALA | HG | `·HG·` | 14 | ALA not a mercury residue |
| ALA | HG1 | `·HG1` | 14 | |
| CYS | HG | `·HG·` | 14 | |
| HG | HG | `HG··` | 13 | A1 |
| MMC | HG | `HG··` | 13 | A4 |
| MMC | HG1 | `HG1·` | 13 | A4 |
| CMH | HG | `HG··` | 13 | A4 |
| EMC | HGX | `HGX·` | 13 | A4 |
| MBO | HG | `HG··` | 13 | A4 |
| HGB | HG | `HG··` | 13 | A4 |
| BE7 | HG | `HG··` | 13 | A4 |
| PMB | HG | `HG··` | 13 | A4 |
| PMB | H | `·H··` | 14 | 1 char |
| PMB | HH | `·HH·` | 14 | not `HG` |
| DCE | CL1 | `CL1·` | 13 | A5 |
| DCE | CL | `CL··` | 13 | A5 |
| DCE | C1 | `·C1·` | 14 | |
| ALA | CL1 | `·CL1` | 14 | |
| CPT | CL2 | `CL2·` | 13 | A5 |
| 0QE | CL | `CL··` | 13 | A5 |
| HLT | CL3 | `CL3·` | 13 | A5 |
| VGH | CL | `CL··` | 13 | A5 |
| UCL | CL | `CL··` | 13 | A5 |
| XMM | CLX | `CLX·` | 13 | A5 |
| OEC | CA1 | `CA1·` | 13 | A6 |
| OEX | CA1 | `CA1·` | 13 | A6 |
| OEX | CA2 | `·CA2` | 14 | pair not listed |
| PLL | PD | `PD··` | 13 | A6 |
| ALA | PD | `·PD·` | 14 | |
| ALA | CA1 | `·CA1` | 14 | |
| C14 | C14 | `·C14` | 14 | excluded |
| C15 | C15 | `·C15` | 14 | excluded |
| F9F | F9F | `·F9F` | 14 | excluded |
| OAN | OAN | `·OAN` | 14 | excluded |
| BLM | NI | `·NI·` | 14 | excluded (beats A2) |
| BZG | CO | `·CO·` | 14 | excluded |
| BZG | NI | `·NI·` | 14 | excluded |
| VNL | CO1 | `·CO1` | 14 | excluded |
| VNL | CO2 | `·CO2` | 14 | excluded |
| PF5 | FE1 | `·FE1` | 14 | excluded |
| PF5 | FE2 | `·FE2` | 14 | excluded |
| UNL | UNL | `·UNL` | 14 | excluded |
| XYZ | UNL | `·UNL` | 14 | |
| ALA | UNK | `UNK·` | 13 | A3 |
| UNK | UNK | `UNK·` | 13 | A1/A3 |
| UNL | UNK | `UNK·` | 13 | A3 |
| C14 | C1 | `·C1·` | 14 | |
| C14 | `C14 ` (trailing space) | `C14·` | 13 | 4 chars |
| ALA | (empty) | `····` | blank | empty |
| ALA | X | `·X··` | 14 | |
| HOH | O | `·O··` | 14 | |
| HOH | OW | `·OW·` | 14 | |
| SOL | HW1 | `·HW1` | 14 | |
| ALA | ABCDE | `ABCD` | 13 | truncated to 4 |
| ABCDE | ABCDE | `ABCD` | 13 | truncated |
| ALA | F | `·F··` | 14 | |
| ALA | CAL | `·CAL` | 14 | |
| CAL | CAL | `CAL·` | 13 | A1 |
| ALA | NA+ | `·NA+` | 14 | |
| ALA | `C 1` | `·C·1` | 14 | inner space kept |
| ALA | ` C` (leading space) | `··C·` | 15 | the name is not stripped |
| ALAX | CA | `·CA·` | 14 | |
| CYS | SG | `·SG·` | 14 | |
| MET | SD | `·SD·` | 14 | |
| ALA | K | `·K··` | 14 | |
| K | K | `·K··` | 14 | 1 char beats A1 |
| ALA | IOD | `·IOD` | 14 | |
| IOD | IOD | `IOD·` | 13 | A1 |
| ALA | 1H | `·1H·` | 14 | |
| ALA | 1HA | `·1HA` | 14 | |
| DA | OP1 | `·OP1` | 14 | |
| DA | C1' | `·C1'` | 14 | |
| DA | H5'' | `H5''` | 13 | 4 chars |
| DA | HO5' | `HO5'` | 13 | 4 chars |
