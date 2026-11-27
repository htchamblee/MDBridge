# PDB entries

| File | Origin | Retrieved | Licence |
|---|---|---|---|
| `4jrd_atoms.tsv` | the ATOM records of `_atom_site` (both RNA strands; waters and ammonium ions left out) of https://files.wwpdb.org/pub/pdb/data/structures/divided/mmCIF/jr/4jrd.cif.gz (Last-Modified 2025-09-25 17:17:49 GMT; latest revision 2024-02-28), as the file writes them | 2026-11-26 | CC0 1.0 (wwPDB data usage policy, https://www.wwpdb.org/about/usage-policies; attribution as a courtesy) |
| `4jrd_entry.tsv` | the entry's id, title, method, resolution, deposition date, latest revision and primary citation | 2026-11-26 | CC0 1.0 |
| `triplex_levels.tsv` | the stacked U·A-U and C+·G-C base triples of the four RNA triple helices below (entry, level, triple, purine / Watson-Crick / Hoogsteen residue as author chain and number), chosen as below | 2026-11-27 | facts (a list of residues) |
| `{3p22,4plx,6svs,7jnh}_triplex_atoms.tsv` | the ATOM records of `_atom_site` of the residues of those triples (the fields of `4jrd_atoms.tsv`, as the files write them) of https://files.wwpdb.org/pub/pdb/data/structures/divided/mmCIF/{p2,pl,sv,jn}/{id}.cif.gz | 2026-11-27 | CC0 1.0 |
| `{3p22,4plx,6svs,7jnh}_entry.tsv` | each entry's id, title, method, resolution, deposition date, latest revision and primary citation | 2026-11-27 | CC0 1.0 |

PDB 4JRD: crystal structure of the parallel double-stranded helix of poly(A) RNA, X-ray, 1.00 A;
Safaee, N., Noronha, A.M., Rodionov, D., Kozlov, G., Wilds, C.J., Sheldrick, G.M., Gehring, K.,
*Structure of the parallel duplex of poly(A) RNA: evaluation of a 50 year-old prediction*, Angew. Chem.
Int. Ed. 52 (2013) 10370-10373, doi:10.1002/anie.201303461 (citation key `pdb-4jrd`).

Cut with `python3 -I tools/extract_pdb_entry.py 4jrd.cif.gz`, which checks the file's SHA-256
`fede84f06b2771793f447a598756aaab40c0f35a2611fdf8ede2bb016acbe408` (36 611 bytes). Used by
`tools/gen_chem_tables.py` for the structure builder's 4JRD nucleotide (chain A residue 10) and the
screw of the duplex (docs/DIFFERENCES.md CHD-02).

## RNA triple helices (docs/DIFFERENCES.md CH-18)

The reference geometry of the builder's Hoogsteen duplex and pyrimidine-motif triplex
(`tools/gen_chem_tables.py` `triplex_tables` -> `src/chem/triplex_tables.inc`):

| Entry | Structure | Method, resolution | Primary citation (key) |
|---|---|---|---|
| 3P22 | KSHV PAN ENE with A9 RNA | X-ray, 2.50 A | Mitton-Fry et al., *Poly(A) tail recognition by a viral RNA element through assembly of a triple helix*, Science 330 (2010) 1244-1247, doi:10.1126/science.1195858 |
| 4PLX | the MALAT1 ENE + A-rich tract | X-ray, 3.10 A | Brown et al., *Structural insights into the stabilization of MALAT1 noncoding RNA by a bipartite triple helix*, Nat. Struct. Mol. Biol. 21 (2014) 633-640, doi:10.1038/nsmb.2844 |
| 6SVS | a U·A-U-rich triple helix of 11 consecutive triples | X-ray, 2.50 A | Ruszkowska et al., *Molecular structure of a U·A-U-rich RNA triple helix with 11 consecutive base triples*, Nucleic Acids Res. 48 (2020) 3304-3314, doi:10.1093/nar/gkz1222 |
| 7JNH | a double ENE with an A28 RNA | X-ray, 2.89 A | Torabi et al., *RNA stabilization by a poly(A) tail 3'-end binding pocket and other modes of poly(A)-RNA interaction*, Science 371 (2021) eabe6523, doi:10.1126/science.abe6523 |

The triples of `triplex_levels.tsv` were found by H-bond criteria (own scratch code of the
2026-11-27 geometry check, docs/reports/2026-11-27-non-wc-geometry.md): a purine with its
Watson-Crick partner (A N1-U N3, G N1-C N3 and the second H-bond under 3.7 A) and a pyrimidine on its
Hoogsteen edge (U N3-A N7 and U O4-A N6, C N3-G N7 and C N4-G O6 under 3.7 A); kept are the triples
of runs in which two consecutive triples have their purine, Watson-Crick and Hoogsteen residues linked
O3'-P to the next triple's (the Hoogsteen strand parallel to the purine strand), all of them standard
A / G / C / U residues (the modified 3'-terminal residues of 4PLX and 6SVS left out): 60 U·A-U and 5
C+·G-C triples, 52 steps (model 1, no alternate locations in these residues).
`tools/gen_chem_tables.py` re-checks the triples' H-bonds and finds the steps by those links.

Cut with `python3 -I tools/extract_pdb_triplexes.py 3p22.cif.gz 4plx.cif.gz 6svs.cif.gz 7jnh.cif.gz`,
which checks each file's SHA-256 (the files of 2026-11-27; Last-Modified 2025-09-25 for 3P22, 4PLX and
6SVS, 2026-08-16 for 7JNH):

| File | SHA-256 | Bytes |
|---|---|---|
| 3p22.cif.gz | `22aca08798b54da2619dcc588b7b2211f98c23815a2dcd8b376c1b01a8ff885c` | 186 392 |
| 4plx.cif.gz | `ab735099ead9988d39993742e91e24013e1dd283728974b8ae4934f98946d4f3` | 248 808 |
| 6svs.cif.gz | `57b79d3b391209907832d42412fd120f60674429dd10a9a9ec8664fe26d6a65d` | 187 583 |
| 7jnh.cif.gz | `f0494fb55ac4f815d6e8264e0e2b57f651f35d6d9240edc77feadecac4f31abc` | 79 058 |

The NMR structures of the telomerase RNA pseudoknots (PDB 2M8K, 2K96) were used as checks of the
measured geometry (the report), not as sources.
