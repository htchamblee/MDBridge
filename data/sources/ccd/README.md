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
| `hg_cl_atom_names.tsv` | every released component with an atom named `HG` whose element is Hg, or named `CL` whose element is Cl |
