# microMg / nanoMg Mg²⁺ parameter sets (extracts not committed: licence question open)

MDATK 0.0.4.1 (`mdatk/workflows/mg_params.py`) copies the optimised 12-6 Mg²⁺ sets of Grotz,
Cruz-León & Schwierz "verbatim from the authors' GROMACS parameter files" at three commits of
github.com/bio-phys. Owner decision D5 of `docs/specs/m7-plan.md` is to extract them from those
files into this folder. `tools/extract_mg_params.py DIR` does that (`mg_sets.tsv`: ion, chloride,
Mg²⁺-Cl⁻ pair per (model, water); `mg_rna_pairs.tsv`: the 27 Mg²⁺-RNA pairs per set; sigma nm,
epsilon kJ/mol, as the files print them, with file and line), but **its output is not committed**:
none of the three repositories states a licence.

| Repository | Commit (MDATK's) | Files read (SHA-256) | Licence found |
|---|---|---|---|
| https://github.com/bio-phys/Magnesium-FFs (TIP3P) | `e71379d1a88b54221018af53ca54d3c18066ba05` (2021-12-02, the head on 2026-11-28) | `ff_Mg.itp` `08e7349831324e2c34e6917f5a7eecefaf8c301a7c4056c4218d7cbb05999ddb` | none: no LICENSE / COPYING file at any commit; README: "If you use our optimized parameters for Magnesium please refer to: K. K. Grotz, S. Cruz-Leon, N. Schwierz (JCTC, 2021 ...)", "USE AT YOUR OWN RISK!!"; the file header: "For implementation please cite ..." |
| https://github.com/bio-phys/optimizedMgFFs (SPC/E, TIP3P-FB, TIP4P/2005, TIP4P-Ew, TIP4P-D) | `1df53e9c94fc14a97864ac35c6951085bb60d1c0` (2021-12-09, the head) | `ff_Mg_spce.itp` `7188362e...` , `ff_Mg_tip3pfb.itp` `52a7c973...`, `ff_Mg_tip4pew.itp` `1a88d2f7...`, `ff_Mg_tip4pd.itp` `71cde72b...`, `ff_Mg_tip4p2005.itp` `a8b29ee0...` (full sums in the tool) | none (README: citation request, "USE AT YOUR OWN RISK!!") |
| https://github.com/bio-phys/MgFF_OPC (OPC: the SPC/E values) | `cf8148a3131c198b19f3a2896b68bad6f6f62f5b` (2021-12-09, the head) | `ff_Mg_opc.itp` `8f84f1f7...` | none (README: citation request, "USE AT YOUR OWN RISK!!") |

Retrieved 2026-11-28 (`git clone`; GitHub reports no licence for any of them). Without a licence
the files are all rights reserved, so they are not redistributed here, and the extracts wait for
the owner's decision.

## What was checked (with the files outside the repository)

* MDATK 0.0.4.1's `MG_SETS` equals the files exactly: all 14 sets (micro / nano x TIP3P, SPC/E,
  OPC, TIP3P-FB, TIP4P-Ew, TIP4P-D, TIP4P/2005), ion, chloride, Mg²⁺-Cl⁻ and the 378 Mg²⁺-RNA pairs as
  doubles, the pair order and the source URLs (`tools/gen_workflow_tables.py --mg-extracts DIR`).
* Every Mg²⁺-RNA pair of a set is the Lorentz-Berthelot pair of the ion and the file's own copy of
  the parmBSC0χOL3 type scaled by one λσ and one λε per set, to the digits printed, and these factors
  equal the λσ RNA / λε RNA of Table 1 of the papers in all 14 sets.

## Recording them as facts from the papers instead

The three papers are open access under CC BY 4.0 (Crossref licence records; Europe PMC full text):
Grotz, Cruz-León & Schwierz, J. Chem. Theory Comput. 17, 2530 (2021), doi:10.1021/acs.jctc.0c01281
(PMC8047801); Grotz & Schwierz, J. Chem. Theory Comput. 18, 526 (2022), doi:10.1021/acs.jctc.1c00791
(PMC8757469); Grotz & Schwierz, J. Chem. Phys. 156, 114501 (2022), doi:10.1063/5.0087292. Their
Table 1 gives σ, ε of the ion, the λ factors and σ, ε of the Mg²⁺-Cl⁻ and Mg²⁺-OP pairs, but rounded
(σ 0.1019 nm, ε 235.80 kJ/mol where the file has 0.101934 / 235.804213166223), and with the
transcription errors MDATK's docstring lists (TIP4P-D microMg σ(Mg-OP), TIP3P nanoMg σ(Mg-OP), OPC
microMg σ(Mg-Cl), the TIP4P/2005 nanoMg Mg-OP pair). The λ factors (4 decimals) agree with the files;
the chloride models of Table 2 agree too. So the papers alone cannot give MDATK's numbers (and so
its frcmod and ParmEd lines) byte for byte: the full-precision ion parameters are only in the files.
Options for the owner: (a) treat the numbers as facts and commit the extracts (as
`../forcefields/amber_lj_types.tsv` does for AmberTools' parameter files, citing the papers); (b) ask
the authors for a licence or permission; (c) use the papers' rounded values with the λ derivation
(a deliberate difference from MDATK in every Mg²⁺ file).
