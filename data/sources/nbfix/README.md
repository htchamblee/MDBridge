# NBfix19, 0BPh, CP and TIP4P/2005 parameters

Published force-field numbers the Amber step types write, each with its source. Hand-typed from the
sources named (no file copied); the numbers are facts, the publications are cited when a step uses
them (MDATK's citation table: `steinbrecher2012`, `kuhrova2019`, `mlynsky2023`, `abascal2005`).

| File | Content | Used by |
|---|---|---|
| `nbfix_pairs.tsv` | the NBfix19 and 0BPh pair terms (Rmin_ij in Angstrom, epsilon_ij in kcal/mol, as the ParmEd `changeLJPair` lines write them), for the OL3 atom types | `tools/gen_workflow_tables.py` -> `src/workflow/tables.inc` (`kNbfixPairs`, the ParmEd lines and the L-RNA pairs, checked equal to MDATK 0.0.4.1's `NBFIX19_LINES`, `BPH0_LINES`, `lrna_pair_lines`) |
| `cp_types.tsv` | the CP phosphate oxygen types OP / OR: Rmin/2, epsilon, and the mass / polarizability the cp-lib helper writes for them | `tools/extract_workflow_helpers.py` checks the cp-lib helper's `CP_NONBOND` / `CP_MASS` against it |
| `tip4p2005.tsv` | the TIP4P/2005 water model's constants | `tools/extract_workflow_helpers.py` checks the tip4p2005 helper's constants against it |

| Source | Status | Licence |
|---|---|---|
| 0BPh: Mlýnský, Kührová, Stadlbauer, Krepl, Otyepka, Banáš & Šponer, *Simple Adjustment of Intranucleotide Base-Phosphate Interaction in the OL3 AMBER Force Field Improves RNA Simulations*, J. Chem. Theory Comput. 19, 8423-8433 (2023), doi:10.1021/acs.jctc.3c00990 (Europe PMC PMC10687871, read 2026-11-28): "between H5–OR/OS and H4–OR/OS atom types, by 0.25 Å to 2.8808 and 2.9308 Å ... Depths of the potential well (ε i,j parameters) were kept at a default value of 0.0505 kcal/mol" | checked | CC BY 4.0 (article) |
| NBfix19: MDATK v0.0.4.1 `mdatk/workflows/amber.py:79-87` (the owner's table), which attributes it to "Mlynsky V., Kuehrova P., Bussi G., Otyepka M., Sponer J., and Banas P., 2026" (no journal or DOI) as an NBfix replacement of gHBfix19 (Kührová et al., J. Chem. Theory Comput. 15, 3288 (2019), doi:10.1021/acs.jctc.8b00955) | **publication not found**: the 2019 article (author manuscript PMC7491206) discusses NBfix but gives none of these numbers; Europe PMC and GitHub searches (2026-11-28) found no 2026 article or file with them besides MDATK. Source to be confirmed: the owner expects the numbers in the Supporting Information of Kührová et al. 2019 (`kuhrova2019`) and will supply it; until then the table stays as MDATK has it | to be confirmed (the SI of kuhrova2019) |
| CP: Steinbrecher, Latzer & Case, *Revised AMBER Parameters for Bioorganic Phosphates*, J. Chem. Theory Comput. 8, 4405-4412 (2012), doi:10.1021/ct300613v, Table 1 | cited, not re-read (not open access); the same Rmin/2 and epsilon are in AmberTools 26 `dat/leap/parm/frcmod.phosaa10` (NONB lines 73 OP, 75 OR; the phospho-amino-acid types of the same paper), checked 2026-11-28 (package of `../amber/README.md`). The masses / polarizabilities are those of O2 (0.434) and OS (0.465) in `parm10.dat` MASS lines 52 and 54, which the CP types replace (frcmod.phosaa10 gives OR 0.434; the column has no effect in non-polarizable runs) | facts |
| TIP4P/2005: Abascal & Vega, *A general purpose model for the condensed phases of water: TIP4P/2005*, J. Chem. Phys. 123, 234505 (2005), doi:10.1063/1.2121687, Table I (and the TIP4P geometry of section II) | cited, not re-read here (not open access); the values are those MDATK's helper cites | facts |
