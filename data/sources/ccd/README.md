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
