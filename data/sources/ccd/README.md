# wwPDB Chemical Component Dictionary extracts

- **Origin:** https://files.wwpdb.org/pub/pdb/data/monomers/components.cif.gz
  (file last modified 2026-11-14 17:40:24 GMT), retrieved 2026-11-24.
- **Licence:** CC0 1.0 Universal (public domain dedication); wwPDB data usage policy,
  https://www.wwpdb.org/about/usage-policies. Attribution is given as a courtesy.
- **Cut with:** `python3 -I tools/extract_ccd.py components.cif.gz`.

| File | Content |
|---|---|
| `standard_residue_atoms.tsv` | comp_id, atom_id, alt_atom_id (PDB v2.3 spelling), element, `pdbx_backbone_atom_flag`, `pdbx_leaving_atom_flag` of every atom of the residues in `../pdb_format/standard_residues.tsv` |
| `monatomic_components.tsv` | every released (REL) component that has exactly one atom: comp_id, atom_id, element, formal charge, name |
| `protein_list_residue_atoms.tsv` | atom_id, alt_atom_id, element, name and type of every atom of the released components whose id is a non-standard residue name of `tools/gen_data_tables.py`'s protein list (computed from the committed force-field and PDB-format extracts). Only the peptide-linking ones (HIP, HSE, HYP, LYN) and their atom_id are used, as reference names of the collision rule; the others merely share an id with a force-field residue (ASH, HSD, ... are unrelated chemicals in the CCD) |
| `l_deoxyribonucleotide_atoms.tsv` | comp_id, base (one-letter code), atom_id, standard_atom_id and element of every atom of the released `L-DNA LINKING` components whose one-letter code is A, C, G or T and whose atoms and bonds match those of DA, DC, DG or DT (the mirror images 0DA, 0DC, 0DG, 0DT). standard_atom_id is the matching atom of the standard residue: `tools/extract_ccd.py` matches the two components' bond graphs (elements and bonds kept), keeping names wherever the bonds allow, so it differs only where the CCD names an atom its own way (0DA's N9A ... HN61, H2P; 0DT's C5M) |
| `hg_cl_atom_names.tsv` | every released component with an atom named `HG` whose element is Hg, or named `CL` whose element is Cl |
| `ribonucleotide_components.tsv` | comp_id, name, type, formula, formal charge, one-letter code, modification date and InChI of the standard ribonucleotides A, C, G, U and the L-ribonucleotides 0A, 0C, 0G, 0U |
| `ribonucleotide_atoms.tsv` | every atom of those eight components: element, formal charge, aromatic / leaving flags, CIP label, ordinal, model and ideal coordinates, and standard_atom_id / standard_ordinal, the matching atom of the standard residue (the atom itself for A ... U). The match is `match_atoms`' bond-graph match; the two hydrogens of a CH2 / NH2 group that the CCD names its own way (0A's HN6 / HN6A, H5' / H5'A) are named by their position in the ideal coordinates (`name_by_geometry`): HN6, HN6A, H5', H5'A are the mirror images of H62, H61, H5'', H5'. The L components are checked to be the mirror images (CIP labels, chiral volumes of the ideal coordinates) |
| `ribonucleotide_bonds.tsv` | every bond of those eight components: atom ids, order, aromatic flag, ordinal |

The ribonucleotide extracts are the reference of the restraint dictionaries and the hydrogen naming
(`tools/gen_chem_tables.py`, `src/chem/component_tables.inc`). The CCD does not agree with itself
everywhere on which of two such hydrogens is which: the model coordinates of A, C, G, U put H5' on
the pro-R position, their ideal coordinates (and the IUPAC-IUB 1983 rule, H5' = H5'1 = pro-S, and
the deposited NMR entries 2KOC and 1ZIH) on the pro-S one; 0G's ideal coordinates put H21 where its
own model coordinates, and G's and DG's (model and ideal), put H22. The generator takes the
standard residues' ideal coordinates as the reference and prints the others as notes.
