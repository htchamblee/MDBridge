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

The triples of `triplex_levels.tsv` are chosen by `python3 -I tools/select_pdb_triplexes.py 3p22.cif.gz
4plx.cif.gz 6svs.cif.gz 7jnh.cif.gz` (the files below, SHA-256 checked; standard library only): a
standard purine with its Watson-Crick partner (A N1-U N3 and A N6-U O4, G N1-C N3 and G O6-C N4) and
a pyrimidine on its Hoogsteen edge (A N7-U N3 and A N6-U O4, G N7-C N3 and G O6-C N4), every one of
these donor-acceptor distances under 3.7 A, per purine the triple whose six distances add up to the
least (a stacked neighbour of a partner can come within the cutoff too); kept are the triples that
take part in a step, two triples whose purine, Watson-Crick and Hoogsteen residues are linked O3'-P
(under 2.0 A) to the next triple's (the Watson-Crick strand antiparallel), with every atom of the
nucleosides the tables use (C1', base, sugar ring with O2', O3', C5'), standard A / G / C / U residues
only (model 1, alternate location `.` or `A`; the modified 3'-terminal residues of 4PLX and 6SVS end
their runs): 61 U·A-U and 5 C+·G-C triples, 53 steps. (The table of 2026-11-27 had been chosen by scratch
code that missed 4PLX C75·C38-C15; the tool reproduces it with that triple added.)
`tools/gen_chem_tables.py` re-checks all four H-bonds of every triple (3.7 A) and finds the steps by
those links.

Statistics (`tools/gen_chem_tables.py` `triplex_tables`; docs/DIFFERENCES.md CH-20): the copies of
one molecule in a crystal (the four ENE·A9 complexes of 3P22, the three MALAT1 molecules of 4PLX, the
two of 6SVS) repeat the same levels to 0.1-0.2°; they are counted once. A level's *site* is its entry
with its Watson-Crick and Hoogsteen residue numbers (the copies share them), a step's site is its first
level's: 66 levels at 30 sites (28 U·A-U, 2 C+·G-C), 53 steps at 24 sites. Means are weighted by 1 / copies of
the site; standard deviations are those of the sites' means (screw parameters, base-plane separation,
inclination, H-bonds), the pose spreads the RMS over the sites. Leaving out the end-of-run steps (a
run's first and last step: 12 of the 24 sites remain) moves the triplex twist from 30.50° to 31.05°,
rise 3.021 to 2.997 A, base planes 3.263 to 3.242 A (spreads 1.84°, 0.17 A, 0.07 A); leaving out the
steps with a lattice contact (any heavy atom of a triple within 4.0 A of another molecule of the
crystal, symmetry mates included: 15 sites remain) moves them to 30.57°, 2.985 A, 3.251 A, both
together (7 sites) to 31.35°, 2.974 A, 3.239 A; the duplex (purine and Hoogsteen strands) likewise
(30.57° / 2.979 A / 3.270 A against 31.52° / 2.954 A / 3.268 A, 30.38° / 2.955 A / 3.267 A, 31.13° /
2.939 A / 3.268 A; spreads 3.21°, 0.22 A, 0.08 A). No mean moves by its spread, so every site is used.
(The lattice contacts were found with gemmi 0.7.5 in a scratch check, not part of the generator.)

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
