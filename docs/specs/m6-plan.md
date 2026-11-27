# Milestone 6 plan: the `chem` module

Goal (ARCHITECTURE "Migration order" 6 + "Definition of done"): MDATK v0.0.4.1's chemistry layer in
C++: PDB hygiene (L-RNA detection from the sugar geometry, residue-name schemes, chains, renumbering,
SEQRES / FASTA, 5'-phosphates), base-pair annotation (REMARK 950, dot-bracket), the structure builder
(ideal A-form helices of either hand, the other pairings and triplexes, base mutations, the rigid-body
fit into a map), restraint dictionaries (SMILES, renaming acedrg / eLBOW output to the CCD atom names)
and the density fit the analysis window and `render` need (MDATK's `density fit` / `simulate`, the
*Global search* / *Refine* buttons that milestone 5 shows disabled, M5 D2, and `auto_fit`, which
milestone 4 left out, M4 D9). With them the command-line groups `pdb`, `build`, `dict` and `density
fit` / `simulate`. PDB hygiene stays in C++ and milestone 7's Lua generators call it
(ARCHITECTURE "Modules").

Licensing: MDATK code (owner's) is ported freely. No MDAnalysis behaviour is involved (section 3), so
no clean-room specs are needed. New third-party material (an FFT, SciPy's Powell and least-squares
solvers, helix templates, CCD components) is listed in section 3 with its licence; the owner decisions
of section 7 name the choices.

References: MDATK v0.0.4.1 = `htchamblee/MDATK` tag `v0.0.4.1` (9875bb6); "0.1.1" = its later `main`
(1efa0ef, version 0.1.1), read for sections 2 and 5. Line numbers are v0.0.4.1's unless a version is
named. MDBridge base: c8e07b8 (milestone 5 complete).

---------------------------------------------------------------------------------------------------

## 1. Inventory of MDATK 0.0.4.1's chem layer

Sizes: `pdbtools.py` 1 229 lines, `annotate.py` 231, `secstruct.py` 422, `builder.py` 2 324,
`dictionaries.py` 1 445, `density.py` 718, `gui/pdbdialog.py` 298, `gui/builderdialog.py` 424. None of
them imports MDAnalysis; `builder.py` and `density.py` use NumPy and SciPy (`least_squares`, `cKDTree`,
`fft`, `ndimage`, `optimize.minimize(Powell)`), `annotate.py` NumPy; `pdbtools.py` and
`dictionaries.py` are standard library only (both are also copied into generated scripts, section
4.5).

"MDBridge" columns: **done** = ported with parity tests (milestone and DIFFERENCES id), **part** =
some of it, **-** = missing.

### 1.1 PDB hygiene (`pdbtools.py`)

| Function | MDATK (file:line) | What it does | MDBridge |
|---|---|---|---|
| data model `Atom`, `Residue`, `Structure` | pdbtools.py:89-145 | fixed-column atoms (columns 28-76 kept verbatim as `tail`), residues, TER segments, header / trailer / CONECT / MODEL | done (M2: `io::PdbAtom`, `PdbResidue`, `PdbStructure`, src/io/include/mdbridge/io/pdb.hpp) |
| `parse_residue_name` | :155-173 | A, A5, RA, LA, 0A, DA, ADE ... -> base, sugar, hand, terminal | done (M2; L-DNA names NA-01 follow 0.1.1) |
| `scheme_name` | :176-186 | name of a nucleotide in `amber` / `ccd` / `rosetta` / `pdb` (L: LA LA5 LA3 / 0A 0A5 0A3 / 0A; single L residue LA5) | done (M6 W1-A: `chem::schemeName`, src/chem/pdbtools.cpp) |
| `_nucleotide`, `_is_polymer`, `_is_solvent`, `_code`, `_chirality` | :189-218 | backbone-guarded name recognition, polymer / solvent tests, one-letter code | done (M6 W1-A: `chem::nucleotide`, `isPolymer`, `isSolvent`, `residueCode`, `residueChirality`; CH-05, CH-06) |
| `read_pdb` (+ `_fused_field`, `_tokens_tail`, `_parse_atom`) | :221-381 | ChimeraX repair (`GC:`), first model, ANISOU, SIGATM / SIGUIJ / MASTER dropped, notes | done (M2: `readPdbText`, `readPdbFile`; tests/data/io/pdb/pdb_cases.json) |
| `write_pdb` (+ `_element`, `_atom_name_field`, `_residue_field`, `_id_columns`, `_atom_line`, `_seqres_insert_index`, `_remap_conect`) | :384-483 | columns, serials (TER takes one), TER per segment, END, SEQRES placement, CONECT strip / remap | done (M2: `writePdb`, `seqresLines`) |
| `chiral_volume`, `sugar_chirality` | :486-522 | signed volumes at C4', C3', C1' (`L_CHIRAL_SIGNS`, :61-63), `MIN_CHIRAL_VOLUME` 0.5 Å^3, all usable centres must agree, a flat centre prevents a call | done (M6 W1-A: `chem::chiralVolume`, `sugarChirality`) |
| `parse_ranges`, `_in_ranges` | :525-556 | `C:45-73`, `C:45`, `C`, `:45-73`, lists, tuples | done (M6 W1-A: `chem::parseRanges`, `inRanges`, `RangeSpec`; CH-15) |
| `assign_chirality` | :579-630 | override ranges (forced L) > detection > name; the notes (L sugars in D names, D sugars in L names, forced against geometry, undetermined, L-DNA without names, ranges that match nothing) | done (M6 W1-A: `chem::assignChirality`) |
| `_terminals`, `convert_names` | :633-695 | 5' / 3' / single-residue states per TER segment; renaming with the notes (terminal names off the ends, 5'-phosphates on Amber 5' units, CCD 5' ends that keep the component, the single L residue) | done (M6 W1-A: `chem::segmentTerminals`, `convertNames`; CH-06) |
| `to_lrna` | :698-706 | the `to_lrna.sh` replacement: force L in ranges, rename only those | done (M6 W1-A: `chem::toLrna`) |
| `set_records` | :709-728 | `keep`, `rosetta` (polymer residues ATOM), `pdb` (non-standard names HETATM) | done (M6 W1-A: `chem::setRecords`) |
| `_resegment`, `_polymer_segments`, `restore_chains` | :733-795 | chain ids from a reference PDB (matched by base letter; `copy_numbers`), a list per TER segment, or A, B, C ... | done (M6 W1-A: `chem::restoreChains`, `restoreChainsFromReference`) |
| `split_chains`, `combine_chains`, `order_chains` | :798-894 | split every further segment of a chain or at residues (file order), combine (duplicate-number warning), order (alphabetical or given; solvent last) | done (M6 W1-A: `chem::splitChains`, `combineChains`, `orderChains`; CH-03) |
| `strip_5prime_phosphates` | :897-934 | P/OP1/OP2/OP3 off each segment's first nucleotide, lone phosphate residues removed, HO5' placed 0.96 Å from O5' when the residue has hydrogens | done (M6 W1-A: `chem::strip5primePhosphates`) |
| `renumber_residues` | :937-970 | continuous or per chain, from `start`, `keep_gaps`, insertion codes cleared, out-of-order and > 9999 warnings | done (M6 W1-A: `chem::renumberResidues`; CH-01) |
| `_sequence_chains`, `seqres_records`, `_compress`, `fasta` | :973-1027 | SEQRES (13 names per line) and FASTA (RNA / DNA / L-RNA / protein / mixed, L positions) | done (M6 W1-A: `chem::sequenceChains`, `seqresNames`, `seqresRecords`, `fasta`; CH-07) |
| `CleanResult`, `clean` | :1030-1097 | the whole chain of operations in MDATK's order, stale REMARK 950 dropped, CONECT / 99 999-atom / SEQRES notes, record policy following the scheme | done (M6 W1-A: `chem::clean`, `PdbSource`; CH-04, CH-13) |
| `summary` | :1100-1114 | `pdb check`'s JSON (chains, numbers, sequence, names, chirality counts, other residues, CONECT count, notes) | done (M6 W1-A: `chem::summary`, `chem::check`) |
| `add_commands`, `clean_options`, `run`, `main` | :1117-1229 | the `clean`, `check`, `to-lrna` parsers and runner (also the stand-alone copy's `main`) | done (M6 W2-C: src/app/commands/chem_commands.cpp; the stand-alone copy: M7) |

### 1.2 Base pairs, dot-bracket and annotation (`secstruct.py`, `annotate.py`)

MDBridge already has, from milestone 3: `chem::dotbracketToPairs`, `pairsToDotbracket` (levels,
letters, chain breaks), `stems`, `crosses`, `uniquePairs`, `baseLetter`, `canonicalName`, `lwKind`,
`pairRules`, and the builtin detector `chem::PairDetector` (canonical Watson-Crick / G·U pairs from
H-bond distances, MDATK's `PAIR_RULES`, `SCORE_TIE`, C1' window; src/chem/secstruct.cpp, bit-exact
with `-ffp-contract=off`), used by the `basepairs` analysis (src/analysis/basepairs.cpp), the sequence
panels and the example dataset; tests/data/chem/secstruct_cases.json (make_secstruct_cases.py), NA-01
for L-DNA and the backbone guard. Barnaba (`BarnabaRunner`, secstruct.py:277-410) is GPL and never
used: DIFFERENCES "No Barnaba".

| Function | MDATK | What it does | MDBridge |
|---|---|---|---|
| `annotate` | annotate.py:145-231 | REMARK 950 records: header, `ENGINE`, `PAIRS n (WATSON-CRICK, WOBBLE, OTHER)`, `PAIR` per pair, `STACK` (Barnaba only), `CHAINS` / `SEQ` / `DBN` per strand group wrapped to 60, notes (models, ChimeraX lines copied unchanged, dropped pairs, gaps); the original text kept verbatim | done (M6 W1-B: `chem::annotate`, src/chem/annotate.cpp; CH-14, D13) |
| `_res_field`, `_wrapped`, `_strands`, `_strand_groups`, `_dotbracket_records`, `_gap_notes` | :54-133 | record fields, strands = TER segments of nucleotides, groups of strands joined by pairs | done (M6 W1-B: `chem::annotationResidueField`, `annotationWrapped`, `annotationStrands`, `annotationStrandGroups`, ...) |
| `_insert_records` | :136-142 | REMARK 950 / MASTER removed, records before the first MODEL / ATOM / HETATM | done (M2: `io::insertHeaderRecords`) |

### 1.3 Structure builder (`builder.py`)

| Part | MDATK | What it does | MDBridge |
|---|---|---|---|
| Arnott A-RNA monomer, base templates, hydrogen rules | builder.py:84-135, 137-243 | cylindrical coordinates from AmberTools' `nucgen.dat` `$ARNA` (section 3.3); `_place`, `_hydrogens` (ideal C-H / N-H / O-H, HO2' torsion 180°), `template`, wobble turns G −13.3° / U 11.4°, `_screw` (32.7°, 2.81 Å), `_dyad` | done (M6 W1-D: src/chem/helix.cpp; CHD-01) |
| Sequences and messages | :244-457 | `parse_sequence` (T read as U, `5'-...-3'`), `complement`, `helix_kind`, `check_pairing`, `check_triplex`, `auto_partner`, `resolve_partner`, `pair_types`, `pairing_kinds`, `resolve_third`, `triplex_kinds`, `check_helix` (the levelled messages the planner and the window show) | done (M6 W1-D: `chem::parseSequence`, `complement`, `helixKind`, `checkPairing`, `checkTriplex`, `resolvePartner`, `checkHelix`, ..., src/chem/helix.cpp) |
| Watson-Crick helix | :460-690 | `_tail` (column limits), `_strand_residue` (5'-OH or 5'-phosphate with OP3, HO5' / HO3'), `_parse_pair` / `_parse_ids` (chains, start numbers), `build_helix`: duplex / single strand, D or L (x -> −x mirror), names through `convert_names` (default scheme `rosetta`-like: A C G U / 0A 0C 0G 0U, no terminals), `REMARK 951` header, mismatch / wobble notes, centring | done (M6 W1-D: `chem::buildHelix` with pairing `wc`, names through W1-A's `chem::convertNames` / `setRecords`; CH-01, CH-02, D13) |
| Other pairings and triplexes | :693-1878 | PDB 4JRD nucleotide (:820-829, CC0), `PAIRINGS` (reverse-wc, hoogsteen, reverse-hoogsteen, hh) and triplex motifs, H-bond tables (:850-865), `_pair_fit` (least squares, :1010-1063), the reference helix fit `_helix_fit` (:1363-1375; `least_squares` from `_HELIX_START`, cached per kind), `_level_pose` (mismatches moved off, :1624-1669), `_close_levels` (joint closure: rigid levels + per-nucleotide torsions, sparse-Jacobian TRF with LSMR, contact terms refreshed with `cKDTree.query_pairs`, :1453-1581), contacts and geometry notes (`_close_contacts`, `_geometry_notes`, `_torsion_note`, :1746-1878), `_turn_terminal_hydroxyls` (:693-720) | done (M6 W2-A: src/chem/helix_fit.cpp behind `chem::buildHelix`; CH-16 continuation closure, CH-17 continuous torsion term, NU-03) |
| Mutations | :1881-2117 | `parse_mutations` (`A:5G`, `A:5>G`, `:5G`), `mutated_name` (scheme, hand, terminal kept), `_partners` (`PairDetector`), `_swap_base` (weighted Kabsch on C1' + glycosidic frame, B-factor / record from C1', altlocs merged, chirality re-checked), `mutate` (compensatory partners, `--wobble`, REMARK 950 dropped) | done (M6 W2-B: `chem::parseMutations`, `mutatedName`, `mutate`, src/chem/mutate.cpp; CH-08, CH-09) |
| Rigid-body fit | :2120-2250 | `model_atoms` (heavy atoms, first altloc), `transform_pdb_text` (coordinates and ANISOU only), `fit_to_map` (global / local through `density.fit_model`; damaged columns rewritten through `write_pdb`), notes (CC, rotation, shift, CRYST1 kept) | done (M2: `io::transformPdbText`, `io::newAtomTail`; M6 W2-B: `chem::atomElement`, `modelAtoms`, `fitToMap`, src/chem/fit_to_map.cpp; CH-10) |
| Command line | :2253-2324 | `helix` / `mutate` / `fit` parsers, `run` writes the output | done (M6 W2-C: chem_commands.cpp) |

### 1.4 Restraint dictionaries (`dictionaries.py`)

| Part | MDATK | What it does | MDBridge |
|---|---|---|---|
| Reference data | dictionaries.py:72-131, 135-259 | atoms / bonds / charges / hydrogen names of the CCD's 0A 0C 0G 0U (as in the CCP4 monomer library), MDATK's L SMILES, the eight 5' / 3' variants (`Component`: leaving atoms, `formula`) -> `COMPONENTS` (:1439) | done (M6 W1-C: `chem::dictionaryComponents`, src/chem/dictionaries.cpp) |
| SMILES and graphs | :262-631 | `MolGraph`, `parse_smiles` (brackets, ring closures, `%nn`, aromatic, stereo `@` / `@@`, charges, implicit H), `_all_matches` (deterministic backtracking, `MAX_MATCHES` 2000), `_score`, `match_graph`, stereo parity, `describe_mismatch`, `smiles_mapping`, `check_smiles`, `identify_smiles`, `get_component` | done (M6 W1-C: src/chem/smiles.cpp, dictionaries.cpp) |
| mmCIF | :634-857 | tokens, `parse_cif`, `cif_value`, `write_cif` | done (M2: src/io/cif.cpp; tests/data/io/cif) |
| Renaming | :857-1251 | `finish_dictionary` / `rename_dictionary`: component block, dictionary graph, heavy-atom isomorphism to the reference names, hydrogen naming (H5'/H5'' and H61/H62 from coordinates, `_hpos_cost`), phosphate protonation, chirality check from the dictionary's coordinates and `_chem_comp_chir`, leaving-atom removal (OP3 of the 3' variants), `_chem_comp` id / name / group, CCD-style input left unchanged | done (M6 W1-C: `chem::finishDictionary`; DIC-01, DIC-02) |
| Files | :1254-1290 | `is_dictionary_file`, `find_dictionary` | done (M2: `io::isDictionaryFile`, `io::findDictionary`) |
| PubChem | :1293-1339 | PUG REST lookup (`SMILES`, then `IsomericSMILES`), only on request | done (M6 W2-C: src/app/commands/pubchem.cpp, Qt Network transport src/app/pubchem_qt.cpp; D7) |
| Command line | :1342-1436 | `list [--json]`, `smiles CODE`, `finish INPUT... --code --output [--group --name --source-id --allow-protonation --keep-names]`, `pubchem COMPOUND [--code --timeout]` | done (M6 W2-C) |

### 1.5 Density fit and simulated maps (`density.py`, `pipeline.py`)

| Part | MDATK | What it does | MDBridge |
|---|---|---|---|
| `RigidTransform`, `anchored`, `anchor_points`, `rotvec_matrix`, `quaternion_matrix`, `rotation_angle` | density.py:59-170 | placements and their algebra | done (M4: `core::RigidTransform`, src/core/rigid_transform.cpp; `render::anchored`, RN-32) |
| `rotation_samples` | :147-162 | super-Fibonacci rotations (Alexa 2022), at least 24, step ≥ 3° | done (W1-E: `chem::rotationSamples`, src/chem/density_fit.cpp) |
| `simulate_map`, `splat_gaussians`, `atom_weights`, `load_map_cached`, `map_level` | :47-56, 176-259 | molmap-style Gaussians (sigma factor 1/(π√2)) | done (M2: `io::simulateMap`, `atomWeights`, `loadMapCached`, `mapLevel`) |
| `sample_map`, `pearson` | :249-273 | trilinear sampling (`map_coordinates` order 1, constant NaN), correlation about the mean | part (order-1 sampling in src/render/molecule/ndimage.cpp, mode "nearest") |
| `FitProblem` | :285-470 | envelope points (5 % of the maximum, ≤ 60 000), `cc` with the inside fraction (full score ≥ 0.5, zero ≤ 0.25), `at_resolution` (Gaussian low-pass copy), `flips` (180° about the principal axes), `refine` (Powell, 6 parameters, `direc` = step·I, xtol / ftol / maxfev), `global_search` (rotations × FFT translations, two peaks per rotation, distinct poses (12, ≤ 3 per site), coarse-to-fine at 8 Å, best 3 + flips of the best 2 refined, message) | done (W1-E: `chem::FitProblem`, several threads, `core::Reporter`) |
| `_CoarseSearch`, `_rfftn`, `_irfftn`, `_fft_size` | :480-590 | coarse padded grid (≤ ~96 points per axis, 2^a 3^b 5^c sizes), anti-alias smoothing, Roseman's locally normalised correlation by FFT, peak suppression | done (W1-E: `chem::CoarseSearch`, `chem::fftSize`) |
| `density_centroid`, `centre_transform` | :593-608 | the map's centre above a level; *Centre on model* | done (`numerics::densityCentroid`, W0-C, used by the window's *Centre on model*; M6 W1-E: `chem::centreTransform`) |
| `fit_model` | :611-634 | global search, or refinement from a start (`manual` / `identity`) | done (W1-E: `chem::fitModel`, `chem::centreTransform`) |
| shadow (`project_density`, `shadow_rgba`, ...) | :637-718 | the density shadow | done (M4: src/render/molecule/density_shadow.cpp, RN-30..33) |
| `Pipeline.density_frame_atoms` | pipeline.py:365-390 | superposed atoms of the fit frame (the `atoms` and `trace` analyses on `ref` and the frame) and the reference anchor | done (M6 W2-C: `app::densityFrameAtoms`, src/app/commands/density_pipeline.cpp) |
| `Pipeline.fit_density` | pipeline.py:392-436 | start (given, stored + anchor, or centred on the density), fit, store `transform`, `anchor`, `fit_cc`, `fit_info` | done (M6 W2-C: `app::fitDensity`) |
| `Pipeline.ensure_density_fit` | pipeline.py:438-455 | fit before rendering when a shown map has no placement (`auto_fit`); a failed fit is a warning | done (M6 W2-C: `app::ensureDensityFit` in `render` / `preview`; M4 D9's warning retired) |
| `DensityParams`, `density_problems` | config.py:280-305, 578-600 | the `[density]` settings and their checks | done (M2 / M5: `core::Project`, src/core/project.cpp:633) |

### 1.6 Command line: MDATK commands vs MDBridge

MDBridge offers none of these yet (src/app/commands/commands.hpp:6: "density fit / simulate, pdb,
build, dict, workflow follow with later milestones"), except `density info` (M2). Every command below
takes `--project` / `-p` / `--params` (the project whose step log records the step; `dict`: `finish`
only) and logs `step` / `title` / `command` / `files` / `summary` / `citations` with `logged()`
(cli.py:494-571); notes go to stderr prefixed `pdbtools:`, `annotate:`, `build:`, `dictionaries:`;
errors print `error: <message>` with exit status 1 (cli.py:998-1020), usage errors exit 2 (argparse).

| Command (cli.py) | Options (MDATK) | Step log | MDBridge |
|---|---|---|---|
| `pdb clean INPUT` (892-897; pdbtools.py:1117-1150) | `-o/--output` (default stdout), `--scheme {amber,ccd,rosetta,pdb}`, `--records {keep,rosetta,pdb}`, `--no-terminal-names`, `--keep-conect`, `--no-detect`, `--l-ranges`, `--chain-reference`, `--chain-ids`, `--assign-chains`, `--copy-numbers`, `--split-segments`, `--split-at`, `--combine`, `--combine-into`, `--sort-chains`, `--renumber`, `--start` (1), `--per-chain`, `--keep-gaps`, `--seqres`, `--strip-5p`, `--fasta FILE` | `pdb clean`, "PDB clean-up of NAME" | done (M6 W2-C) |
| `pdb check INPUT` | `--l-ranges` | none; JSON on stdout (`summary` + notes) | done (M6 W2-C) |
| `pdb to-lrna INPUT RANGES` | `--scheme` (amber), `--no-terminal-names`, `-o` | `pdb to-lrna` | done (M6 W2-C) |
| `pdb annotate INPUT` (898-909) | `-o` (required), `--no-pairs`, `--no-stacking`, `--no-dotbracket`, `--engine {auto,builtin,barnaba}` | `pdb annotate`, "Base-pair annotation of NAME", citation `barnaba` when used | done (M6 W2-C) |
| `build helix` (builder.py:2259-2287) | `--seq/--sequence` (required), `--partner` (auto), `--chirality {D,L}`, `--pairing {wc,reverse-wc,hoogsteen,reverse-hoogsteen,hh}`, `--triplex {none,hoogsteen,reverse-hoogsteen}`, `--third-strand/--third` (auto), `--chains` (A,B), `--start` (1,1), `--scheme`, `--no-terminal-names`, `--phosphate5`, `--hydrogens`, `-o` (required) | `build helix`, "Ideal A-form D-RNA helix SEQ" (`helix_kind`), citations `arnott-arna` (+ `pdb-4jrd` for other pairings) | done (M6 W2-C; the other pairings and triplexes: W2-A) |
| `build mutate INPUT MUTATIONS` | `--no-compensatory`, `--wobble`, `-o` | `build mutate`, "Mutation of NAME: MUTATIONS" | done (M6 W2-C) |
| `build fit INPUT MAP` | `--resolution` (required), `--search {global,local}`, `--step` (30), `-o`; prints `{"output", "cc"}` | `build fit`, "Rigid-body fit of NAME into MAP" | done (M6 W2-C) |
| `dict list` / `smiles` / `finish` / `pubchem` (cli.py:551-566) | section 1.4 | `finish` only: `dict finish`, "Restraint dictionary CODE", citation `ccd` | done (M6 W2-C) |
| `density fit` (cli.py:436-471, 870-879) | the analysis commands' project options (`_common`), `--map`, `--resolution`, `--level`, `--frame`, `--search {global,manual}`, `--step`, `--save [PRESET]`; JSON `{cc, start, info, transform[, saved]}` | none | done (M6 W2-C) |
| `density simulate` (cli.py:474-491, 880-886) | `_common`, `--resolution` (required), `--spacing`, `--frame`, `--placed`, `-O/--output` (required); JSON `{map, frame, atoms, grid_xyz, voxel_A, frame_of_reference}` | none | done (M6 W2-C) |
| `density info MAP` | | | done (M2) |
| `workflow restraints` (cli.py:983-992) | helix restraints for a parm7 | | milestone 7 (workflow) |

### 1.7 GUI

| Element | MDATK | Milestone |
|---|---|---|
| Structure view page, density section: *Global search*, *Refine* (`fit_density(search)`, a job that runs `Pipeline.fit_density` and stores the placement, then reloads the preview), status line | gui/app.py:924-1056 | **6** (shown disabled since M5, "comes with milestone 6": src/app/gui/pages/structure_page.cpp:20, 129-136) |
| *Centre on model*, *Clear*, map mode (move / rotate the map with the mouse) | app.py:964-1031, canvas.py | done (M5 W2-A / W2-B) |
| *PDB clean-up…* window (`PDBToolsWindow`: input / output, names and chirality, chains, numbering, sequence, annotation, *Check input*, *Write cleaned PDB*, report) | gui/pdbdialog.py:39-298; menu app.py:163, 1632-1638 | D1 (section 7); done (M6 W3-B: src/app/gui/tools/pdb_tools_window.cpp behind tools/host.hpp; DIFFERENCES TW-xx) |
| *Structure builder…* window (`BuilderWindow`: *Ideal helix*, *Mutate*, *Fit to map* pages, *Check*, *Build*, background jobs with progress and cancel) | gui/builderdialog.py:39-424; menu app.py:165, 1640-1646 | D1; done (M6 W3-B: tools/builder_window.cpp; the other pairings with W2-A, TW-07) |

MDATK 0.0.11 moved both windows into RNAdo behind a small host interface (0.1.1 `gui/host.py`); M5 D1
put RNAdo into milestone 7.

### 1.8 MDATK tests that specify the chem layer

| File | Tests | Counterpart |
|---|---|---|
| test_pdbtools.py | 32 (chiral volumes, detection / override / fallback, ranges, names and schemes, round trips, `to_lrna`, ChimeraX repair, columns, models / ANISOU, record policy, CONECT, chain restore / split / combine / order, renumbering, SEQRES, FASTA, `clean` end to end, CLI + step log, embedded copy, terminal names, split order, single L residue, `--strip-5p`) | W1-A (GoogleTest + generated cases); the embedded-copy test in M7 |
| test_annotate.py | 8 (builtin records and placement, options and L-RNA, 60-column wrap, Barnaba records, CLI + step log, real Barnaba, strands from TER segments / ions / verbatim models, clean drops stale annotation / MASTER / missing Barnaba) | W1-B; the Barnaba tests as the DIFFERENCES "No Barnaba" behaviour |
| test_builder.py | 31 (duplex geometry, L mirror image, single strand / validation, partners with wobble and mismatch, names / terminals / hydrogens / phosphate, C3'-endo sugars, backbone vs an experimental A-form helix, sugar stereocentres (RDKit, optional), mutations, L mutations, mutation errors, global and local fits, ANISOU rotation, CLI + step log, planner step (M7), builder window (D1), negative start / empty partner / DNA partner / altlocs, L CCD names vs `pdb clean`, hh vs 4JRD, pairings, triplex, L pairings mirrored, pairing messages, pairing CLI and planner step, no close contacts, reverse-WC bases clear their sugars, C2' vs next C5', closed strands keep the pairs, contact report with hydrogens, tleap (optional), builder window pairings) | W1-D, W2-A, W2-B; planner and window tests with M7 / D1 |
| test_dictionaries.py | 25 (component table, built-in SMILES vs graph and stereo, wrong stereo / compounds, SMILES parser, deterministic matching, RDKit (optional), CIF (done M2), rename round trip, group / name / source id / keep names, wrong compound / D stereo, phosphate protonation, L hydrogen names, CCD-style input unchanged, `find_dictionary` (done M2), PubChem new / old property names and errors, CLI, CLI + step log; 7 planner-step tests (M7)) | W1-C (W2-C for the CLI) |
| test_density.py | 15, of them for M6: rotations cover SO(3), global search recovers the pose, pipeline fit + video with density, `density` CLI, global search ignores unrelated density and validates steps, placement follows the molecule and failures degrade | W1-E, W2-C |
| test_gui.py | `test_density_map_fit_shadow_and_drag` (fit part), `test_pdb_tools_window_checks_and_cleans`, `test_builder_non_wc_helix_builds_in_the_background` | W3-A; W3-B (D1) |
| test_secstruct.py, test_basepairs.py | 9 + the analysis tests | done (M3) |

### 1.9 Behaviour confirmed by running MDATK

Probes (scratch scripts against v0.0.4.1 and 0.1.1, Python 3.13, numpy 2.5.5, scipy 1.18.1):
`build helix --seq GGCAUGCC --chirality L` writes 0G/0C/0A/0U ATOM records with the `REMARK 951`
header; `pdb clean --scheme amber --renumber` renames to LG5 ... LC3 with the notes of section 1.1;
`--scheme ccd --seqres --fasta` writes 0G5 / 0C3, SEQRES after the header and `L-RNA` FASTA headers;
`build mutate ... A:3G` changes the partner (compensatory) and drops the old base's hydrogens;
`pdb annotate` (builtin) writes `PAIRS 8 (WATSON-CRICK 8, ...)`, `SEQ GGGAUGCC&GGCAUCCC`,
`DBN ((((((((&))))))))`; `dict list` lists the 12 codes; `dict pubchem` reaches PubChem. Timings:
a Watson-Crick helix well under a second; the other pairings and triplexes of 8 levels 3-7 s (the
least-squares closures); a global density fit of an 8-bp duplex at 6 Å 3.5 s, of the demo system
4.9 s (551 rotations at 30°). The CH items of section 5 come from these probes.

---------------------------------------------------------------------------------------------------

## 2. MDATK after 0.0.4.1 (to 0.1.1): chem-related changes

`density.py` and `dictionaries.py` are unchanged in 0.1.1 (byte-identical files); `annotate.py`
only gained the DSSR engine; `secstruct.py` the backbone guard and the CCD atom names (MDBridge
NA-01). `pdbtools.py` grew to 2 183 lines, `builder.py` to 2 578; new chem-side modules:
`chimera.py` (2 389), `modelcheck.py` (1 837), `interfaces.py` (1 822), `mapscore.py` (1 290),
`ensemble.py` (2 146), `ensemble_stats.py` (389), `rosetta.py` (1 479), `frames.py` (817), `dssr.py`
(317). Classes: **M6** = adopt now (a fix of a v0.0.4.1 behaviour: a CH item, with the version named
in DIFFERENCES); **later** = with milestone 7 (the generators that need it); **post** = the separate
post-parity phase. Parity tests stay on v0.0.4.1 (M5 D18).

| Change | Version (CHANGELOG / file) | Class |
|---|---|---|
| DNA written with RNA names (ADE, CYT, GUA; pre-2007 A, C, G) is DNA when the sugar is 2'-deoxy | 0.1.0 (pdbtools.py:317-328) | M6 (CH-05) |
| A nucleotide bonded to a chain-end residue the schemes do not know (5'-GTP, modified nucleotide) keeps its internal name | 0.1.0 (pdbtools.py:1364-1390, the O3'-P bond test) | M6 (CH-06; the LINK part post) |
| Mutation specs without a chain read the whole number (`15G` = blank chain, residue 15) | 0.1.0 (builder.py:2094) | M6 (CH-08) |
| SEQRES names have three characters: longer nucleotide names as their CCD component, others UNK with a WARNING | 0.0.17 (pdbtools.py:1861-1895) | M6 (CH-07) |
| L-DNA names read; backbone guard; CCD atom names of 0DA / 0DT read as standard names | 0.0.13-0.0.17 | done (NA-01) |
| Every scheme names DNA (`scheme_name(..., sugar="dna")`, amber DA5 / DA3 / DAN, L-DNA 0DA ... / LDA ...), `pdb clean` notes the DNA it leaves, `write_pdb` writes 0DA / 0DT component atom names, `--dna` | 0.0.13 / 0.0.17 | post (with B-DNA; changes v0.0.4.1 output of DNA files only by notes) |
| Three-character residue names only (L-DNA LDA at chain ends) | 0.0.17 | post (with the DNA schemes) |
| Scheme `standard`, record policy `atom`, `--strip-h`, `--drop-virtual` (Rosetta VRT / XXX), `infer_base` | 0.0.12 | post (the Rosetta decoy and model-check workflows) |
| `pdb check --expect-l / --expect-d / --strict / --table`, `residue_hands` | 0.0.12 | post |
| LINK records follow their residues through every operation; LINK-bonded ends keep internal names; `--strip-5p` keeps a LINKed phosphate; phosphate-only residues out of sequences | 0.0.13 (`follow_links`, pdbtools.py:755-1035) | post (with chimeras) |
| Empty-slot files (`REMARK 951 empty slot`), `read_pdb` refuses them | 0.0.12 | post (model check) |
| B-form DNA helices (`build helix --na dna`; Arnott `$ABDNA`, phosphate turned 6°, G·T wobble), DNA mutations | 0.0.13 | post |
| D/L chimeras (`build chimera --design`, `--link-partners`) | 0.0.13 / 0.0.18 (chimera.py) | post |
| Model check `pdb qc` (clashes, O3'-P links with strain, loop compaction, hands, kinks, gates) | 0.0.12 (modelcheck.py) | post |
| `pdb interface`, `pdb hybrids`, `decoys census / clusters / support`, `rosetta scores / audit`, `frame check / reframe`, `ensemble census / compare / regions`, `map score` (MDATK CC or ChimeraX Q-score) | 0.0.12 (interfaces.py, mapscore.py, ensemble.py, rosetta.py, frames.py) | post (separate domain: decoy ensembles) |
| DSSR annotation engine (`--engine dssr`, x3dna-dssr run as an external program, licensed separately) | 0.0.12 (dssr.py) | post (an external program, never bundled) |
| PDB clean-up and builder windows behind a host interface; *Check* reports a file that is not UTF-8; *Cancel* of a fit logged as cancelled | 0.0.11 (gui/host.py), 0.1.0, 0.0.13 | with the windows (D1) |
| Amber build: L-RNA names checked before tleap, generated L-RNA / L-DNA libraries, embedded PDB tools copy changes | 0.0.8-0.0.19 | later (M7 generators) / post |

`density fit`, `simulate`, `build fit` and the `dict` commands keep v0.0.4.1's options in 0.1.1;
`build helix` adds `--na`; `pdb clean` adds `--strip-h`, `--drop-virtual`, `--dna`; `pdb check` adds
four options; `pdb annotate` adds the `dssr` engine; `pdb qc`, `pdb interface`, `pdb hybrids`,
`build chimera` and the groups `decoys`, `rosetta`, `frame`, `ensemble`, `map` are new (D8).

---------------------------------------------------------------------------------------------------

## 3. Clean-room boundary, data and licences

### 3.1 MDAnalysis

No chem module imports MDAnalysis (section 1). The one path that touches it is the density fit's
frame (`Pipeline.density_frame_atoms`, pipeline.py:365-390): the `atoms` and `trace` analyses on the
reference and fit frames, which MDBridge has since milestone 3 (written from the m3 clean-room specs).
`spatial.kabsch` (mutations) is MDATK's own (NumPy SVD). So milestone 6 needs **no new clean-room
specs**; spec writers are not involved. Implementers may read all of MDATK's chem code and run it, as
for milestones 4 and 5, and never open MDAnalysis.

### 3.2 Data tables (CLAUDE.md: from `data/sources/` only)

| Table | MDATK's source | M6 |
|---|---|---|
| Amino-acid, variant, water, two-letter element, standard-name lists of `pdbtools` (pdbtools.py:65-79) | MDATK's own short lists | ported as in milestone 2 (`src/io/pdb.cpp` already has the element list); owner's code, not MDAnalysis |
| CCD components 0A 0C 0G 0U: atoms, elements, charges, bonds, hydrogen names (dictionaries.py:72-126) | "the CCD entries ... as distributed in the CCP4 monomer library (acedrg 243)" | extracted from the wwPDB CCD (`components.cif`, CC0) by `tools/extract_ccd.py` into `data/sources/ccd/l_ribonucleotide_components.tsv`; the generator checks that it equals MDATK's table (D6). The monomer library (its own licence) is not used |
| L SMILES of the 12 codes | MDATK's own (dictionaries.py:94-126) | owner's code; checked against the CCD graph by the existing `check_smiles` logic |
| Arnott A-RNA fibre monomer (r, φ, z of 12 backbone and 8-11 base atoms) | AmberTools `dat/leap/parm/nucgen.dat`, `$ARNA` (builder.py:84-104) | see 3.3 (D5) |
| PDB 4JRD nucleotide and screw (builder.py:820-829) | PDB entry 4JRD (CC0) | extracted from the wwPDB mmCIF of 4JRD by a tool into `data/sources/pdb/4jrd_*.tsv`; the generator checks MDATK's cylindrical numbers against it |
| Backbone geometry (`_GEOM`: Parkinson et al. 1996 / 4JRD means), torsion targets, H-bond pairs, contact limits, `_HELIX_START` | MDATK's own constants with published values | owner's code; Parkinson et al. cited in the citations data (already has `arnott1972`, `pdb-4jrd`) |
| Atomic numbers of the fit (`ATOMIC_NUMBER`, density.py:38-40) | MDATK's own | done (M2 `io::atomWeights`) |

### 3.3 Licences of algorithms and external data

| Item | Licence / status | Plan |
|---|---|---|
| Arnott A-RNA coordinates via `nucgen.dat` | the numbers are published experimental fibre coordinates (Arnott, Hukins & Dover 1972; Chandrasekaran & Arnott 1989); the file ships with AmberTools under AmberTools' licences (GPL / LGPL / other, per component), so it is not copied | record them as facts, as `data/sources/ATTRIBUTION.md` already does for parm10.dat numbers ("numbers with file and line; the files are not copied"): `data/sources/helix/arnott_arna.tsv` with the publication, `nucgen.dat` file and line cited, checked against MDATK's embedded table (D5) |
| PDB 4JRD | wwPDB data CC0 | extract (3.2) |
| wwPDB CCD | CC0 | extract (3.2) |
| Super-Fibonacci rotations (Alexa, CVPR 2022), Roseman's FLCF (2003), ChimeraX `molmap` / `fitmap` conventions | published algorithms; MDATK's code is the owner's | port MDATK's code; no third-party code |
| FFT (`scipy.fft` = pocketfft C++) | BSD-3 (pocketfft, M. Reinecke) | vendor pocketfft's header-only C++ (`pocketfft_hdronly.hpp`, the code `scipy.fft` is built on; the version SciPy 1.18.1 bundles) with its notice, or MDBridge's own radix-2/3/5 FFT (D3) |
| Powell's method and Brent's line search (`scipy.optimize._optimize`) | BSD-3 | port with notice (as the ndimage port, LICENSING.md) |
| `scipy.optimize.least_squares` (TRF, LSMR, finite-difference Jacobian with sparsity) | BSD-3 (LSMR from Fong & Saunders' reference code, BSD) | port with notice, or another solver (D4) |
| `scipy.spatial.cKDTree` (`query_pairs`, `query_ball_point`) | BSD-3 | MDBridge's cell grid (`analysis/spatial`) with an explicit pair order (risk 4) |
| `scipy.ndimage.gaussian_filter` (per-axis sigma), `map_coordinates` order 1 mode constant | BSD-3 | already ported (src/render/molecule/ndimage.cpp, M4 notice); move to `core` (or a small `numerics` library) so `chem` can use it without depending on `render` |
| Kabsch (NumPy SVD) | MDATK's own; LAPACK results machine dependent in the last bits | MDBridge's existing superposition code (`core::linalg`, tests/data/rigid) |
| Leontis-Westhof classes, the builtin pair rules, `WCc` / `GUc` codes | nomenclature (facts) and MDATK's own rules | done (M3) |
| Barnaba | GPL v3 | never used (DIFFERENCES "No Barnaba") |
| x3dna-dssr (0.1.x engine) | licensed separately, not open source | not in M6; if adopted later, run as an external program only |
| PubChem PUG REST | public web service (PubChem's usage policy: request-rate limits) | on explicit request only, as MDATK (D7) |
| acedrg / eLBOW output | user's files | only parsed |

---------------------------------------------------------------------------------------------------

## 4. Parity strategy

### 4.1 What is compared

| Artifact | Same means |
|---|---|
| PDB text from `clean`, `to-lrna`, `annotate`, `mutate`, a Watson-Crick `build helix` (all options), `transform_pdb_text` | byte for byte (pure Python arithmetic or 3 × 3 NumPy products; where the machine-dependent BLAS / LAPACK path of a 3 × 3 product or of Kabsch flips a third decimal, the regenerate job compares the ATOM records as numbers within 0.0011 Å, as `MACHINE_DEPENDENT` does for the rigid cases) |
| Notes, JSON of `pdb check`, `dict list --json`, `build fit` / `density fit` / `density simulate` JSON | equal text (MDATK's wording, MDBridge's where a CH item changes it) |
| Dot-bracket, SEQ, DBN, PAIR records | equal text (the builtin detector is already bit-exact) |
| Other pairings and triplexes (least-squares fits) | heavy-atom RMSD to MDATK's model ≤ 0.01 Å, every atom ≤ 0.05 Å, names / records / header equal; numbers in notes (twist, rise, fit residual, closest contacts, torsion means) within their printed rounding ± 1 in the last digit (D10) |
| Dictionaries (`finish`) | mmCIF text byte for byte (graph matching is deterministic; hydrogen names from coordinates are discrete choices) |
| Density fit | transform within 0.5° and 0.2 Å of MDATK's, CC within 1e-3, `fit_info` equal after masking the CC digits that differ; `density simulate` maps byte for byte (M2's `simulateMap` already is) |
| Project files after `density fit --save` and GUI fits | MDATK's file after the CH-12 patch, `transform` / `anchor` / `fit_cc` / `fit_info` compared as numbers with the fit tolerance |
| Step-log entries of every command | equal entries (time masked), titles and citations included |

### 4.2 Generators (`tests/parity/`, `regenerate_all.py`)

- `make_pdbtools_cases.py`: MDATK's test fixtures (`tests/test_pdbtools.py` helpers, as
  make_pdb_cases.py loads them) and crafted structures through every function of section 1.1 alone
  (`sugar_chirality` with volumes as `float.hex`, `parse_ranges` texts and errors, `scheme_name` of
  every base × hand × terminal × scheme, `convert_names` per scheme with and without terminals and
  ranges, chain operations, renumbering variants, SEQRES, FASTA, `strip_5prime_phosphates`, `summary`)
  and `clean` over a grid of option combinations (fixed seed), with notes and errors.
- `make_annotate_cases.py`: `annotate` of the fixtures, built helices (D, L, single strand,
  mismatches, several strand groups, ions, models, ChimeraX lines, long sequences for the wrap),
  options, the barnaba engine as MDBridge answers it (CH-14).
- `make_builder_cases.py`: `check_helix` messages over sequences / partners / pairings / triplexes;
  Watson-Crick helices over sequence, partner, chirality, chains, start, scheme, terminals,
  phosphate5, hydrogens (text + sha256); the non-WC builds of MDATK's tests and a grid (text, plus the
  coordinates for the tolerance comparison); mutations (fixtures, L names, altlocs, errors);
  `transform_pdb_text` (done in pdb_cases).
- `make_dict_cases.py`: `COMPONENTS` rows, `parse_smiles` (MDATK's test strings and a fuzzed set),
  `match_graph` / `smiles_mapping` / `check_smiles` / `identify_smiles` results and errors,
  `finish_dictionary` of MDATK's synthetic acedrg / eLBOW dictionaries (test_dictionaries.py builds
  them; make_cif_cases.py already loads that file) with every option.
- `make_density_fit_cases.py`: `rotation_samples` (counts, first / last matrices), `_fft_size`,
  `FitProblem` envelopes and `cc` values, `_CoarseSearch.peaks` for given rotations, `refine` and
  `global_search` / `fit_model` results on simulated maps (helix, demo frame, map with unrelated
  density, partial overlap, the GUI test's map), with NaN-voxel maps patched (CH-11).
- `improvements_chem.py`: the test-side patches of MDATK for the CH items (as improvements_*.py),
  each naming the MDATK version whose design it follows.
- `cli_parity.py`: cases for every command and option of section 1.6 (MDATK's tests' command lines
  first, then one case per option, errors, `-p` step-log entries), `density fit` / `simulate` on the
  demo project, `render` of a project with a map and `auto_fit` (the log's `Density fit: ...` line
  instead of M4's D9 warning). Usage errors: exit status only (`exit_only`), as today.
- `make_gui_cases.py` (M5 W0-A): the fit part of `test_density_map_fit_shadow_and_drag` (*Global
  search* on a simulated map: saved `density` table), and with D1 the two windows' action scripts.

### 4.3 Tests in C++

GoogleTest suites per area (`tests/chem/test_pdbtools.cpp`, `test_annotate.cpp`, `test_builder.cpp`,
`test_dictionaries.cpp`, `test_density_fit.cpp`) read the JSON cases; each MDATK test of section 1.8
gets a counterpart (named after it). Geometry tests MDATK runs with RDKit or tleap become
MDBridge-side checks of the same property (signed chiral volumes at C1', C2', C3', C4'; no close
contacts) without those programs. Qt Test for the window parts (W3).

### 4.4 The API milestone 7's Lua generators call

MDATK's generators call the chem layer while checking and writing plans (0.0.4.1):

| Caller | Functions |
|---|---|
| workflows/amber.py:538-549, 849-858 | `parse_ranges`, `read_pdb`, `assign_chirality` (L ranges checked at plan time); the embedded copy of `pdbtools.py` (4.5) |
| workflows/buildstep.py:134-140, 206, 229-251, 270-305, 311 | `resolve_partner`, `resolve_third`, `pair_types`, `check_helix`, `parse_mutations`, `read_pdb`, `helix_kind`; the scripts run `python -m mdatk build helix / mutate / fit` |
| workflows/dictstep.py:80-86, 133-298, 331 | `read_pdb`, `assign_chirality`, `convert_names(..., "ccd")`, `COMPONENTS`, `DEFAULT_GROUP`, `check_smiles`, `parse_smiles`, `fetch_pubchem_smiles` (only when asked); the embedded copy of `dictionaries.py` |
| workflows/refine.py:137-189, 729 | the REMARK 950 PAIR / STACK fields, `read_pdb`, `_nucleotide`, `assign_chirality`, `secstruct.lw_kind` |

So `chem` exposes these as value-in / value-out functions (strings, lists, small records; notes as
string lists; errors as `io::PdbError` / `BuildError` / `DictionaryError` with MDATK's messages) in
one header, `chem/generator_api.hpp`, which milestone 7 binds to Lua without wrapping logic. M6
tests that header against the generator cases; no Lua in M6.

### 4.5 Embedded Python copies (boundary to milestone 7)

The Amber step copies `pdbtools.py` into the plan folder and runs it with `python3` (amber.py:849-858,
`workflows/scripts/embed`); the dictionaries step copies `dictionaries.py` (dictstep.py:331). Their
text is MDATK's source, so "generated scripts byte for byte" means milestone 7 must either write the
same text (MDATK's code, the owner's) or change the scripts to call `mdbridge pdb clean` / `mdbridge
dict finish` (a deliberate difference). M6's CLI keeps MDATK's options exactly so that either works;
the choice is milestone 7's.

---------------------------------------------------------------------------------------------------

## 5. MDATK chem bugs and fragile behaviour to fix (ids CH-xx)

Found by reading v0.0.4.1 and running it (section 1.9); "0.1.1" says whether MDATK's later `main`
still has it (the same probe run there) and how it was fixed; where it was, MDBridge follows that
design. The bug report of 2026-11-25 (docs/reports/2026-11-25-mdatk-bug-report.md) covers the
analysis, GUI, reader and rendering layers and has no chem items; CH-12 is its MDATK-BUG-006 on the
command line, CH-11 its MDATK-BUG-060 in the fit. Each lands as a DIFFERENCES entry (prefix CH).

| Id | MDATK v0.0.4.1 behaviour (verified) | MDBridge | 0.1.1 |
|---|---|---|---|
| CH-01 | Residue numbers below −999 do not fit columns 23-26: `_id_columns` (pdbtools.py:407-410) writes `A-1002`, shifting every later column, and MDATK's own `read_pdb` then fails (`line 5: unreadable coordinates`). Reached by `build helix --start=-1002,1`, `pdb clean --renumber --start -2000`, `--copy-numbers` | `PdbError` before anything is written: "residue number -1002 (A) does not fit the PDB columns (-999 to 9999)"; numbers above 9999 keep MDATK's wrap (a convention) with its warning | not fixed |
| CH-02 | `build helix` (and `mutate` of such files) wraps numbers above 9999 silently (`--start 9998,1`: `B9999` then `B   0`), while `renumber_residues` warns (pdbtools.py:968-969) | the same WARNING in the build notes | not fixed |
| CH-03 | `--combine-into XY` is accepted; the notes say "into chain XY" but the file has chain X (pdbtools.py:848-869, `[:1]` at :410) | `PdbError("chain ids must be single characters")`, as `restore_chains` (:786) | not fixed |
| CH-04 | `--keep-conect` with repeated input serials (concatenated files, wrapped serials) drops every CONECT record (pdbtools.py:466-468) while the notes say "Kept CONECT records (remapped to the new atom serials)" (:1083-1084) | the note says what happened: "WARNING: dropped N CONECT records: the input's atom serials repeat, so they cannot be remapped" | not fixed |
| CH-05 | DNA written with RNA names (CHARMM ADE / CYT / GUA, pre-2007 A / C / G; C2' and no O2') is renamed as RNA (`G5`, `A3`; pdbtools.py:189-194) | a 2'-deoxy sugar makes it DNA; DNA keeps its names | fixed in 0.1.0 (pdbtools.py:317-328): adopt |
| CH-06 | A segment that starts (ends) with a residue the schemes do not know (PSU, a 5'-GTP) gives the next nucleotide a terminal name (`C5`) and the 5'-phosphate warning (pdbtools.py:633-641) | a nucleotide bonded through O3'-P to that neighbour is internal | fixed in 0.1.0 (`_terminals`, pdbtools.py:1364-1390; without its LINK part): adopt |
| CH-07 | SEQRES writes names longer than three characters in four (`LDA5   C   A   U`; pdbtools.py:985-993), misaligning the record | 0.0.17's rule: a longer nucleotide name as its CCD component, any other as UNK with a WARNING, with a note | fixed in 0.0.17 (`_seqres_names`, pdbtools.py:1861-1895): adopt |
| CH-08 | A mutation without a chain id reads its first digit as the chain: `15G` = chain 1, residue 5 (builder.py:1881) | `15G` = blank chain, residue 15 | fixed in 0.1.0 (builder.py:2094): adopt |
| CH-09 | Mutating a model with hydrogens drops the old base's hydrogens and adds a heavy-atom base (builder.py:1968-2032): a mixed model that tleap and the model checks then treat inconsistently | the new base gets its hydrogens from the builder's ideal rules (`_hydrogens`) when the old base had any, with a note (D11) | not fixed |
| CH-10 | `--resolution nan`, `--spacing -1` (`build fit`, `density fit`, `density simulate`) end in `negative dimensions are not allowed` (density.py:183-189; builder.py:2199 lets NaN through; cli.py:485 reads `--spacing 0` as the default) | "the resolution must be a positive number of Å" / "the grid spacing must be positive" before any work (`--spacing 0` keeps MDATK's default) | not fixed (density.py unchanged) |
| CH-11 | Non-finite voxels in a fit: the Gaussian smoothing of the coarse search and of the coarse-to-fine copy spreads one NaN voxel over 9^3 voxels at a sigma of one voxel (density.py:345, 496-497) and sampling reads them as outside the map; one NaN voxel turns a 0.0° fit into 0.1° and changes the CC | non-finite voxels carry no density (RN-30 / RN-33's rule), a note counts them | not fixed |
| CH-12 | `density fit --save` (cli.py:451-471) writes every path of the project absolute (topology, trajectories, the map), as the GUI save does (AP-06, MDATK-BUG-006) | saved through the edit model (ED-02): paths as the file had them; a new `--map` path stored relative to the project folder when it was given relative | not fixed (verified) |
| CH-13 | A missing input file is reported as `no ATOM/HETATM records found (not a PDB file, or the path does not exist)`: a name that is not a file is parsed as PDB text (pdbtools.py:221-228) | the command line names the missing file (MDBridge's usual file error); text input to the API unchanged | not fixed |
| CH-14 | `pdb annotate` without Barnaba says "install barnaba (pip install barnaba)" (annotate.py:204-205); `--engine barnaba` fails with that hint (annotate.py:184-185, secstruct.py:277-284) | `auto` is the builtin engine; `barnaba` runs the builtin engine with the note `WARNING: Barnaba is not available in MDBridge; using the builtin detector` (printed `annotate: WARNING: ...`, as the analysis does, DIFFERENCES "No Barnaba"); the stacking note without the pip hint | n/a (MDBridge has no Barnaba) |
| CH-15 | A reversed range (`A:4-1`) matches nothing, reported only as "the L ranges match no nucleotide" (pdbtools.py:525-549) | `PdbError("residue range A:4-1 ends before it starts")` | not fixed |

Not bugs, kept: the default scheme of `build helix` (A C G U / 0A 0C 0G 0U, ATOM, no terminals),
`pdb clean` printing to stdout without `-o`, the note texts and prefixes, `REMARK 950` / `951` numbers
and layout, 4-character residue names in columns 18-21, the wrap of numbers above 9999 and serials
above 99 999 (with the warnings), the 2-strand / 3-strand `--chains` defaults, `MAX_MATCHES` 2 000,
the fit's defaults (30°, 12 poses, 60 000 points, 8 Å coarse stage), the local refinement's helical
register ambiguity (a 3.7 Å offset of an 8-bp duplex refines to the screw-related optimum, CC 0.786:
a property of the problem; `fit_to_map`'s docstring notes the related two-fold ambiguity of
self-complementary duplexes).

---------------------------------------------------------------------------------------------------

## 6. Work breakdown

Code layout: `src/chem/` gains `pdbtools.{hpp,cpp}` (on `io::PdbStructure`), `annotate.*`,
`helix.*` (templates, Watson-Crick builder), `helix_fit.*` (other pairings, closures), `mutate.*`,
`smiles.*`, `dictionaries.*`, `density_fit.*`, `generator_api.hpp`; numerics shared with `render`
(`ndimage`, FFT, Powell, least squares) go to `src/core/numerics/` (or a `mdbridge_numerics`
library); `src/app/commands/` gains `pdb.cpp`, `build.cpp`, `dict.cpp` and the `density fit` /
`simulate` commands next to `density info`. `chem` links `core` and `io` only.

### Wave 0 (independent; start together)

| Package | Content | Acceptance |
|---|---|---|
| **W0-A Generators** | section 4.2 generators and `improvements_chem.py` (CH patches), `cli_parity.py` cases for pdb / build / dict / density fit / simulate (skipped by `needs` until the commands exist), regenerate-job steps | `regenerate_all.py --only <steps>` byte-identical twice; MACHINE_DEPENDENT handling for the fitted outputs documented |
| **W0-B Data sources** | `tools/extract_ccd.py` gains the 0A 0C 0G 0U components; a tool for PDB 4JRD (mmCIF, CC0); `data/sources/helix/arnott_arna.tsv` as facts (D5); `tools/gen_chem_tables.py` -> `src/chem/*_tables.inc`; ATTRIBUTION rows; the generator checks MDATK's embedded tables equal them | tables regenerate byte-identically; equality with MDATK's numbers asserted in the generator |
| **W0-C Numerics** | move `ndimage` and `densityCentroid` to core (per-axis Gaussian, order-1 `map_coordinates` mode constant with cval); FFT (D3: pocketfft vendored or own) `rfftn` / `irfftn` 3-D; Powell + Brent port (SciPy 1.18.1, BSD notice); `least_squares` TRF with LSMR, `x_scale`, sparse finite-difference Jacobian (D4); a pair search with cKDTree's `query_pairs` order or a documented sorted order | each against SciPy on recorded cases: FFT ≤ 1e-12 relative (bit-exact if pocketfft), Powell / TRF iterates within 1e-9 on smooth test problems, pairs equal as sets |

### Wave 1

| Package | Needs | Content | Acceptance |
|---|---|---|---|
| **W1-A PDB hygiene** | (W0-A cases) | section 1.1: chirality, ranges, schemes, `convert_names`, `to_lrna`, records, chains, 5'-phosphates, renumbering, SEQRES / FASTA, `clean`, `summary`; CH-01, CH-03..07, CH-13, CH-15 | `make_pdbtools_cases` all equal; the 32 test_pdbtools counterparts (embedded copy: M7) |
| **W1-B Annotation** | W1-A (`_nucleotide`), M3 detector | section 1.2, CH-14 | `make_annotate_cases` equal; 8 counterparts |
| **W1-C Dictionaries** | W0-B | SMILES parser, graphs, matching, `finish` / `rename`, `check` / `identify`, `COMPONENTS` from the CCD extract | `make_dict_cases` equal; 16 counterparts (CIF and `find_dictionary`: M2; the 7 planner-step tests: M7) |
| **W1-D Watson-Crick helix** | W0-B, W1-A (`convert_names`, `set_records`) | templates, hydrogens, wobble, `check_helix` and the sequence functions, `build_helix` for `wc`, header and notes; CH-01 / CH-02 for the builder | WC cases byte for byte; geometry tests (C1'-C1' 10.9 Å, O3'-P 1.59 Å, C3'-endo, mirror image) |
| **W1-E Density fit** | W0-C | `rotation_samples`, `FitProblem`, `_CoarseSearch`, `refine`, `global_search`, `fit_model`, `centre_transform` reuse; CH-10, CH-11; the same result for any thread count (sums in a fixed order, as the M3 runner); cancellation and progress through `core::Reporter` | `make_density_fit_cases` within 4.1's tolerance; test_density fit counterparts; ≤ MDATK's time on the demo |

### Wave 2

| Package | Needs | Content | Acceptance |
|---|---|---|---|
| **W2-A Other pairings and triplexes** | W0-C (TRF, pairs), W1-D | 4JRD nucleotide, `_pair_fit`, `_helix_fit` (reference helices; the refined starts recorded by the generator), `_level_pose`, `_close_levels`, contacts, torsion and geometry notes, `_turn_terminal_hydroxyls` | non-WC cases within D10's tolerance; the pairing and triplex counterparts of test_builder (`test_hh_duplex_matches_4jrd` to `test_contact_report_includes_hydrogens_and_close_heavy_atoms`, 11 tests; the planner part, tleap and window tests excepted); time ≤ MDATK's |
| **W2-B Mutations and fit-to-map** | W1-A, W1-D, W1-E | `parse_mutations` (CH-08), `_swap_base`, compensatory partners, CH-09 (D11); `model_atoms`, `fit_to_map` with the damaged-column rewrite | mutation cases byte for byte (Kabsch within 0.0011 Å); fit cases within tolerance |
| **W2-C Commands and pipeline** | W1-A..E, W2-B | `pdb clean / check / to-lrna / annotate`, `build helix / mutate / fit`, `dict list / smiles / finish / pubchem` (D7), `density fit` (with `--save` through the edit model, CH-12) / `simulate`; step-log entries and citations; `Pipeline` fit frame, `fit_density`, `ensure_density_fit` in `render` / `preview` (M4 D9's warning goes); `generator_api.hpp` | every new `cli_parity.py` case green; MDATK's CLI tests' counterparts |

### Wave 3

| Package | Needs | Content | Acceptance |
|---|---|---|---|
| **W3-A Density fit in the window** | W2-C, M5 job runner | *Global search* and *Refine* enabled (the shell text goes), jobs with progress and cancel, the placement stored as one undo step, status text, preview reload | the fit part of `test_density_map_fit_shadow_and_drag` (saved `density` table within tolerance); cancel leaves the project unchanged |
| **W3-B PDB clean-up and builder windows** (D1) | W2-C, W3-A pattern | Qt windows behind a host interface (0.1.1 `gui/host.py` design) opened from the analysis window's menu until M7 moves them to RNAdo; 0.1.x's window fixes (UTF-8 check, cancel logged as cancelled) | `test_pdb_tools_window_checks_and_cleans`, `test_builder_non_wc_helix_builds_in_the_background` counterparts; GUI action scripts equal MDATK's files |
| **W3-C Docs** | all | DIFFERENCES "Chemistry" section with every CH item (what, why, how tested, version followed), README / BUILDING (new dependency if D3 vendors pocketfft), LICENSING rows (pocketfft, SciPy optimize ports), ATTRIBUTION rows, ARCHITECTURE status, `docs/specs/m6-status.md` | docs reviewed against the code |

### Critical path and parallelism

`W0-C → W1-E → W2-B / W2-C → W3-A` (numerics decide the density fit, which the window and `render`'s
`auto_fit` need) and `W0-C → W2-A` (the least-squares closures are the largest and riskiest port;
W1-D can ship Watson-Crick helices before it). W1-A is on every command's path but has no numeric
risk. Three agents in wave 0, five in wave 1, three in wave 2, two or three in wave 3.

Packages needing owner decisions first: W0-B (D5, D6), W0-C (D3, D4), W1-A (D12, D14), W1-B (D9,
D13), W1-D (D13), W1-E (D2), W2-A (D10), W2-B (D11), W2-C (D7, D8), W3-B (D1).

---------------------------------------------------------------------------------------------------

## 7. Owner decisions (one line each, with recommendation)

All decided by the owner on 2026-11-26: every recommendation below is accepted.

- **D1 GUI scope**: recommend M6 delivers the library, the CLI and the density fit buttons of the
  analysis window (W3-A), and builds the PDB clean-up and builder windows (W3-B) behind a host
  interface opened from the analysis window's menu as in v0.0.4.1, for M7 to move into RNAdo (M5 D1);
  if capacity is short, W3-B moves to M7 and the menu shells say "comes with milestone 7".
- **D2 Density fit algorithm and defaults**: recommend a faithful port of MDATK's (super-Fibonacci
  rotations, FFT locally normalised correlation, Powell refinement, coarse-to-fine at 8 Å, flips;
  defaults 30°, 12 poses, 60 000 points) with tolerance parity, deterministic for any thread count.
- **D3 FFT**: recommend vendoring pocketfft's header-only C++ (BSD-3, the code SciPy 1.18.1 runs),
  which gives SciPy's numbers, over an own radix-2/3/5 FFT.
- **D4 Least-squares solver for the other pairings**: recommend porting SciPy's `least_squares` TRF
  with LSMR and the sparse finite-difference Jacobian (BSD-3, with notice), so the closures follow
  MDATK's iterates, over an own or Eigen Levenberg-Marquardt (results would differ beyond the
  tolerance).
- **D5 Helix geometry source**: recommend recording the Arnott A-RNA monomer as published facts in
  `data/sources/helix/` (publication and `nucgen.dat` file and line cited, the file not copied, as for
  parm10.dat) and extracting 4JRD from the PDB (CC0), each checked equal to MDATK's numbers.
- **D6 Dictionary reference data**: recommend extracting the 0A 0C 0G 0U atoms, bonds and hydrogen
  names from the wwPDB CCD (CC0) rather than taking MDATK's table, whose stated source is the CCP4
  monomer library; the generator checks they agree.
- **D7 PubChem lookup**: recommend `dict pubchem` in M6 through Qt Network in the `mdbridge` program
  only (never in `chem`, only on explicit request, tested against a local stand-in server).
- **D8 0.1.x commands now**: recommend adopting only the fixes (CH-05..08) in M6 and leaving every
  new 0.1.x command and option (section 2: `pdb qc`, decoys, rosetta, frame, ensemble, map,
  chimera, B-DNA, `standard`, `--strip-h`, `--drop-virtual`, LINK following, DSSR) to the
  post-parity phase.
- **D9 Annotation engines**: recommend the builtin engine for `auto` and `barnaba` (with the
  warning, CH-14) and no STACK records in M6; STACK records from MDBridge's own stacking analysis
  can come post-parity.
- **D10 Parity tolerance for fitted structures**: recommend the tolerances of section 4.1 (non-WC
  helices ≤ 0.01 Å heavy-atom RMSD, ≤ 0.05 Å per atom; fits ≤ 0.5° / 0.2 Å, CC ≤ 1e-3) and byte
  parity everywhere else.
- **D11 Hydrogens on mutated bases (CH-09)**: recommend placing the new base's hydrogens when the
  old base had any (a difference beyond 0.1.1), with a note.
- **D12 Non-standard residues**: recommend v0.0.4.1's pass-through (unknown residues keep name and
  record) plus CH-06, without 0.1.x's `infer_base`.
- **D13 Program name in written records**: recommend `REMARK 950 MDBridge BASE-PAIR ANNOTATION` and
  `REMARK 951 MDBridge IDEAL ...` (as CMD-08 names MDBridge in cpptraj scripts), every other byte as
  MDATK writes it; record numbers and layout unchanged (milestone 7's refinement reads the PAIR
  records by keyword).
- **D14 mmCIF coordinates**: recommend PDB only in M6 (as MDATK), with CH-01 / CH-02 for the column
  limits; mmCIF model output post-parity.

Decided by the owner on 2026-11-26 after wave 0 (W0-B found MDATK's tables disagree with their
primary sources; docs/DIFFERENCES.md "Chemistry reference data"):

- **D15 Guanine's amino hydrogens (CHD-01)**: follow the wwPDB CCD, H21 trans to N1 (N3 side) and
  H22 cis to N1 (the Watson-Crick edge), in the builder's hydrogen rules (built and mutated G), the
  0G dictionary rule and the pair hydrogen-bond tables (wc and reverse-wc G·C, reverse-Hoogsteen
  G·G use H22). The expected data patches MDATK the same way (improvements_chem.py).
- **D16 The 4JRD helix frame (CHD-02)**: use the 4JRD nucleotide in the frame of the fitted 4JRD
  screw (twist 41.748411°, rise 3.791704 Å; src/chem/helix_tables.inc), not MDATK's tilted, rounded
  frame, for every non-Watson-Crick build; the expected data patches MDATK with
  `gen_chem_tables.jrd_nucleotide()`, and the hh helix's REMARK 951 prints TWIST 41.7 (MDATK 41.8).

---------------------------------------------------------------------------------------------------


- **D17 Non-Watson-Crick closure (CH-16)** (decided 2026-11-27): MDATK's `_close_levels` (four truncated
  passes from O3'-P gaps up to 6.8 Å, contact pairs refreshed between passes, inexact subproblems) lands
  on different backbone conformers (up to 1.56 Å heavy-atom RMSD, torsions up to 174° apart) under
  1e-14 changes of its start or another BLAS kernel, so D10 against it is unattainable. MDBridge replaces
  it with a continuation closure (from the uniform helix to the level poses in small steps, each solved
  to convergence with TRF's exact subproblem and a sorted, stable contact refresh), which is reproducible
  within D10 and gives better geometry. Owner: "Reliable geometry is infinitely more important than
  runtime"; runtime is not a constraint. The expected data come from MDATK patched the same way.
- **D18 Realistic non-Watson-Crick geometry (CH-18)** (decided 2026-11-27): measured against
  experimental structures (`docs/reports/2026-11-27-non-wc-geometry.md`), only hh's reference helix
  (4JRD) is realistic. MDATK's (0.0.4.1, unchanged in 0.1.1) hoogsteen, reverse-hoogsteen,
  reverse-wc, th and trh helices are unstacked: levels 4.7-8 Å apart against about 3.3 Å.
  - Owner: fix this before W2-A merges.
  - The hoogsteen duplex and the th (pyrimidine-motif) triplex take stacked reference geometry from
    RNA triplex structures: viral and cellular elements such as the KSHV PAN ENE (3P22) and the MALAT1
    ENE (4PLX), with 6SVS and 7JNH, and the telomerase pseudoknots 2M8K and 2K96 as checks.
  - RNA is not known to form the reverse-wc duplex, the reverse-hoogsteen duplex or the trh
    (purine-motif) triplex. RNA builds of these are refused, with a message saying so. They come back
    for DNA and DNA-RNA hybrids (the owner's example: HBV) with DNA structure building.
  - The CH-17 refit of the reference helices is dropped: its reverse-hoogsteen strands interpenetrate
    from 5 levels on. CH-17's smooth term stays in the closure.

## 8. Risks

1. **Least-squares closures (W2-A)**: TRF with LSMR over sparse Jacobians and contact sets refreshed
   between passes is the largest numeric port; small differences in pair order or BLAS sums can send
   a closure to a slightly different optimum. Mitigation: D4 (same algorithm), the pair order of
   cKDTree reproduced or a sorted order applied on both sides (`improvements_chem.py`), D10's
   tolerance, and Watson-Crick helices (most users) shipping first from W1-D.
2. **Density fit reproducibility**: Powell's path depends on every CC value; with near-degenerate
   poses (self-complementary duplexes, helical register) MDATK and MDBridge can pick different
   equivalent poses. Mitigation: parity cases with unique optima, plus symmetry-aware comparisons for
   the documented degenerate cases.
3. **Machine-dependent NumPy paths in "byte-exact" outputs**: 3 × 3 products and Kabsch go through
   BLAS / LAPACK; a last-bit difference can flip a rounded PDB coordinate. Mitigation: the regenerate
   job's numeric fallback (4.1), as for the rigid and render cases.
4. **cKDTree order**: `query_pairs` order feeds the closures' residual order and the contact notes'
   order. Mitigation: W0-C records it; where notes list contacts, check whether MDATK sorts them.
5. **Data provenance** (D5, D6): the Arnott table comes through an AmberTools file; the dictionary
   table through the CCP4 monomer library. Mitigation: primary sources (publications, CCD, PDB) in
   `data/sources/`, legal review notes in ATTRIBUTION.
6. **Scope pressure from 0.1.x**: 0.1.1's chem side is about three times v0.0.4.1's (some 20 000
   lines against 6 400: chimeras, model check, ensembles, map scores). Held by D8; each item's class
   is in section 2.
7. **Embedded Python copies (4.5)**: milestone 7 must decide how generated scripts carry the PDB
   tools and dictionaries code; M6 keeps the CLI compatible with either answer.
8. **GUI windows (D1)**: building them in M6 and re-hosting them in M7 risks doing the layout twice.
   Mitigation: the host interface of 0.1.1 from the start.
9. **Performance**: MDATK's non-WC builds take 3-7 s for 8 levels and its global fit 3.5-5 s; the
   C++ versions should not be slower, and the global fit must stay cancellable from the window.
