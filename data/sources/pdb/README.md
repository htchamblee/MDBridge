# PDB entries

| File | Origin | Retrieved | Licence |
|---|---|---|---|
| `4jrd_atoms.tsv` | the ATOM records of `_atom_site` (both RNA strands; waters and ammonium ions left out) of https://files.wwpdb.org/pub/pdb/data/structures/divided/mmCIF/jr/4jrd.cif.gz (Last-Modified 2025-09-25 17:17:49 GMT; latest revision 2024-02-28), as the file writes them | 2026-11-26 | CC0 1.0 (wwPDB data usage policy, https://www.wwpdb.org/about/usage-policies; attribution as a courtesy) |
| `4jrd_entry.tsv` | the entry's id, title, method, resolution, deposition date, latest revision and primary citation | 2026-11-26 | CC0 1.0 |

PDB 4JRD: crystal structure of the parallel double-stranded helix of poly(A) RNA, X-ray, 1.00 A;
Safaee, N., Noronha, A.M., Rodionov, D., Kozlov, G., Wilds, C.J., Sheldrick, G.M., Gehring, K.,
*Structure of the parallel duplex of poly(A) RNA: evaluation of a 50 year-old prediction*, Angew. Chem.
Int. Ed. 52 (2013) 10370-10373, doi:10.1002/anie.201303461 (citation key `pdb-4jrd`).

Cut with `python3 -I tools/extract_pdb_entry.py 4jrd.cif.gz`, which checks the file's SHA-256
`fede84f06b2771793f447a598756aaab40c0f35a2611fdf8ede2bb016acbe408` (36 611 bytes). Used by
`tools/gen_chem_tables.py` for the structure builder's 4JRD nucleotide (chain A residue 10) and the
screw of the duplex (docs/DIFFERENCES.md CHD-02).
