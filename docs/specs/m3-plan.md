# Milestone 3 plan: the `analysis` module

Goal (ARCHITECTURE "Migration order" 3 + "Definition of done"): every analysis of MDATK v0.0.4.1
(`mdatk/analysis/*.py`) in C++, same numbers within floating-point tolerance on the test
trajectories, `mdbridge analyze` and `mdbridge list` with MDATK's options and output, each MDATK
test of the analyses with a C++ counterpart. The analysis cache is *not* required to be compatible
(ARCHITECTURE "Compatibility"), so cache files are MDBridge's own; keys are kept identical anyway
(cheap, and they pin `params_for` / `selection_texts` / `referenced_groups`).

Licensing (LICENSING.md): MDATK code ported freely; numpy / scipy (BSD) ported with notices
(`THIRD_PARTY_NOTICES.md`, `packaging/licenses/`, LICENSING table); MDAnalysis (LGPL) only through
behaviour-only specs written by a separate spec agent; implementers of those parts never read
MDAnalysis. Barnaba (GPL) is never used.

---------------------------------------------------------------------------------------------------

## 1. Inventory

Registry order (`ANALYSES` dict = import order; used by `list`, by `required_analyses` sorting and
by `expected_keys`): **contacts, energy, hbonds, mg, rmsd, stacking, trace, proximity, atoms,
basepairs**. trace / proximity / atoms / basepairs are `internal` (not user-enabled; pulled in by
panels).

### 1.1 Framework — `analysis/base.py`, `analysis/__init__.py`, `analysis/common.py`

| Piece | Behaviour to port |
|---|---|
| `Analysis` class attributes | `name`, `label`, `version` (default 1; trace 4, atoms 2), `internal`, `group_slots` (0/1/2), `OUTPUT_LABELS` |
| classmethods | `enabled(project)`, `params_for(project)` (dict -> cache key), `selection_texts`, `expected_outputs` ({"series","matrices","pairmaps"}), `groups` (group_a/group_b when 2 slots), `with_groups` (sets group_a/b; RMSD overrides: targets `Heavy` = `group:G and not name H*`, `Backbone` = `group:G and name C4'`), `output_label` (template `{a}`/`{b}` else `"<label>: <out>"`) |
| `cache_key` | `stable_hash({"analysis", "version", "params", "groups": referenced_groups(selection_texts), "inputs": fingerprint, "use_pbc"})` (16 hex). `params_for` for contacts/hbonds/mg/stacking = `dict(vars(section), enabled=None)` (stable_hash drops None -> enabled removed); energy also pops `target_value`; rmsd `{"targets": [vars(t)], reference_frame, superposition, mass_weighted}`; basepairs `{"groups": structure_groups, "ions": mg.mg_selection, engine, stride:int, cutoff:float, wobble:bool, noncanonical:bool}`; trace `{groups, trace_atoms, align, show_ions, ions ("" unless show_ions), group_order}`; proximity `{groups, set, radius, ions}`; atoms `{groups, trace_atoms, align, hydrogens, ions, group_order}` |
| instance ids | `parse_instance("contacts@A+B") -> ("contacts", ("A","B"))`, `instance_id`, `split_key("x@A+B/out")`, `analysis_class` (KeyError `unknown analysis 'x' (known: a, b, ...)`), `instance_project` (ValueErrors: `"<cls> does not take residue groups ('inst')"`, `"'inst': <cls> takes N group(s)"`, `"'inst': unknown residue group(s) A, B"`; returns a *copy* with groups substituted), `key_label`, `expected_keys(project, include_disabled)` |
| `referenced_groups(texts, project)` | regex `group:([A-Za-z0-9_.+\-']+)` + bare group names, transitive through group selections, `[{name, residues, selection}]` sorted by name |
| `AnalysisContext` | `select(text, allow_empty)` (MDATK select with project groups), `group_atoms(name, heavy_only)`, `box(ts)` = `ts.dimensions` only when `analyses.use_pbc`, `warn(msg)` -> log + reporter `raw_log("WARNING: msg")` |
| `run_analyses(ctx, analyses, frame_indices, reporter)` | `log("Preparing <label> ...")` per analysis then `setup`; `log("Analysing N frame(s): <labels joined ', '>")`; iterate `traj[frame_indices]`; per frame `check_cancel`, `times[k] = ts.time`, `process(k, ts)` for each analysis in order; progress `((k+1)/n, "Analysing frame k+1/n")` when `k % max(1, n//200) == 0` or last; returns `([finalize()...], times float64)` |
| `common.py` | `residue_rows(atoms)` (np.unique resindices, searchsorted rows), `without(atoms, other)` (order kept), `labels` = `residue_labels`, `per_row_counts(offsets, rows, values, n_rows)` -> float32 rows x frames via `np.add.at` |
| `AnalysisResult.new_result` | `params = params_for(project)`, `meta = dict(meta, label=label, analysis=name)` (already `core::AnalysisResult`) |

### 1.2 Analyses (settings: `[analyses.*]`, already in `core::AnalysisParams`)

All "per frame" arrays have length n = analysed frames. CSR pair maps (`<p>_offsets` int64 n+1,
`<p>_i`/`<p>_j` int32, `<p>_v` float32) come from `CsrBuilder` (exists). Labels via
`residueLabels` (exists).

**rmsd** (`group_slots=1`). Settings `targets[{label, selection, fit_selection}]`,
`reference_frame`, `superposition`, `mass_weighted`. enabled = `enabled and targets`. Output key per
target `_out_key(label, k)` = `re.sub(r"[^A-Za-z0-9_.+\-]+", "_", label).strip("_") or "target<k>"`.
Setup: error `RMSD: no targets defined`; `ref_pos = clamp(reference_frame, 0, n-1)`; per target
`mob = select(sel, allow_empty)` (warn `RMSD target 'L': 'sel' matched no atoms; its series will
be empty (NaN)`), `fit = select(fit_sel)` if non-blank; weights = masses (float64) when
mass_weighted (`wf` = fit masses, else w); seeks to `frame_indices[ref_pos]` *during setup* and
copies reference positions (float64). Process: NaN stays for empty mob; no superposition ->
`diff = x - ref`; no fit selection -> **`MDAnalysis.analysis.rms.rmsd(x, ref, weights=w,
center=True, superposition=True)`**; fit selection -> `kabsch(xf, ref_fit, wf)`, `diff = (x - mc)
@ R.T + rc - ref`; value = `sqrt(np.average(sq, weights=w))` or `sqrt(sq.mean())`. Finalize: meta
`reference_frame_index=ref_pos`; arrays `rmsd_<n>` float64; series `<out>` label `"RMSD <label>"`
(+ `" (G)"` for an `@G` instance), unit `Å`.

**contacts** (2 slots). `group_a, group_b, cutoff, heavy_only`. Setup: `a, b = group atoms` (heavy
if heavy_only); `same` = equal atom sets; else `b = without(b, a)` (error `contacts: group 'B' lies
entirely inside 'A'`). Process: `pairs_within(a.pos, b.pos, cutoff, box)` on **float32
positions**; if same: keep `a.ix[i] < b.ix[j]` and `|resindex diff| > 1`; `atom_pairs[k] = n
pairs`; residue keys `row_a[i]*n_b + row_b[j]` -> `np.unique(return_counts)` -> `residue_pairs[k]`
and CSR (i, j, count). Finalize: `q_native` (native = keys of first analysed frame;
`np.isin` hits per frame / len(native); zeros if none), `per_res_a/b` = `per_row_counts`.
Arrays: map_*, `atom_pairs` int64, `residue_pairs` int64, `q_native` float64, `per_res_a/b`
float32. Series atom_pairs / residue_pairs / q_native; matrices per_residue_a/b (rows labels_a /
labels_b); pairmap map. Meta labels_a, labels_b, group_a, group_b, cutoff.

**hbonds** (2 slots). `distance_cutoff, use_angle, angle_cutoff, donor_elements, acceptor_elements`.
Setup: `b = without(group b, a)` (error as contacts, "H-bonds: ..."); `ab = a | b` (sorted union);
`side` 0/1; elements `atom_elements(ab)`; residue rows within own side; donors/acceptors by
element (upper-cased csv); bonds from `ab.bonds` (exception or empty -> unknown): H neighbours
`h_of[x]`, heavy-neighbour counts; acceptor N excluded when it carries H or has >= 3 heavy
neighbours. `use_angle = use_angle and has_h` (warn `H-bonds: topology has no hydrogens; using
donor-acceptor distance only`). No bonds but H present: hydrogens assigned on the first processed
frame (`pairs_within(H, cand_donor, 1.25, box)`, nearest donor, strict `<` so the first pair in
pair order wins ties). Donors = candidates with H (when topology has H). Process (positions
float64): for s in 0,1 donors of side s vs acceptors of 1-s within dcut; with angle: expand per
H, `v1 = mi(D-H)`, `v2 = mi(A-H)`, `cos = dot / max(|v1||v2|, 1e-12)`, keep `cos <= cos(radians(
angle_cutoff))`; count; residue pairs (a-row, b-row) unique+counts -> CSR. Arrays map_*, `count`
int64, per_res_a/b float32; meta + `used_angle`, distance_cutoff, angle_cutoff.

**mg** (2 slots). `mg_selection, ligand_selection, ligand_elements, group_a, group_b, inner_cutoff,
outer_cutoff, bridge_shell, geometry, angle_tolerance, occupancy_selection, occupancy_elements`.
Setup: `mg = select(mg_sel, allow_empty)` (warn `Mg: selection 'x' matched no ions; Mg series will
be zero`); `lig = element_filter(without(select(lig_sel), mg), csv)`; group masks (only for known
group names); occupancy residues from `without(select(occ_sel or lig_sel), mg)`, occupancy atoms
element-filtered; `cos90 = cos(radians(90 - tol))`, `cos180 = cos(radians(180 - tol))`. Process
(float64): `pairs_within(mg, lig, outer, box)`; inner = `d <= inner_cutoff`; contacts / outer
counts; bridge shell; **stable argsort by ion**, per ion all ligand pairs `triu_indices`, optional
octahedral test `|c| <= cos90 or c <= cos180`, bridges, cross-group bridges, bridging ions;
occupancy: `pairs_within(mg, occ_atoms, outer)`, per-residue min distance (`np.minimum.at`) ->
class 0/1/2 int8. Arrays contacts, outer, bridges, inter, bridging_ions, inner_residues (int64),
occupancy int8 [n_occ, n]. Matrix `occupancy` kind "class", classes `["none","outer sphere","inner
sphere"]`. Meta occ_labels, n_ions, group_a, group_b.

**stacking** (2 slots, groups only for scope "inter"). `selection, scope, group_a, group_b,
distance_cutoff, angle_cutoff, offset_cutoff, planarity_cutoff`. Rings per residue of `sel.residues`
(`{name: ix}`, last duplicate wins): purine 9 atoms, pyrimidine 6, else `AROMATIC_RESIDUES[resname
.upper()]` (PHE/TYR/TRP/HIS/HIE/HID/HIP/HSD/HSE). Warn `stacking: no aromatic rings found in the
selection`. Side by residue membership. Process (<2 rings: empty CSR): per ring size batch
**`np.linalg.svd`** of centred ring coords -> centroid, normal `vt[-1]`, planarity `s[-1]/sqrt(k)`;
pairs of centroids within distance (inter: rows A vs cols B; all: `i<j` and different residue);
`cos=|ni.nj|`, offsets from projections of `mi(c_j - c_i)`, criteria; count / inter / intra; CSR
(row_of, col_of, 1). Arrays map_*, count, inter, intra int64; series count, inter, intra (unless
scope "inter"); meta row_labels, col_labels, n_rings, scope.

**energy** (2 slots). `cutoff, dielectric, distance_dependent, use_parm7, exclude_bonded`
(`target_value` display only, re-applied by pipeline). Setup: `b = without(...)` (error);
charges when the topology has them (else warn `energy: topology has no partial charges;
electrostatics will be zero`, charge_source "none"); LJ from `LJTable.from_parm7(topology_path)`
when use_parm7 and parm7 and atom count matches (else warn `energy: could not read LJ parameters
from parm7 (<exc>); using generic element parameters`), else `GENERIC_LJ` by element (default
(1.9, 0.1)); `_bonded_exclusions` BFS (depth = exclude_bonded) -> sorted keys `a_local*n_b +
b_local`. Process: pairs within cutoff, drop excluded, `d = max(d, 0.5)`, Coulomb 332.0522,
eps or eps*r, LJ A/r^12 - B/r^6 (generic: rmin = ra+rb, eps = sqrt(ea eb), A = eps rmin^12,
B = 2 eps rmin^6), sums; net force `|sum(rvec * fmag)|`. Arrays total, elec, vdw, force float64,
n_pairs int64; series total (target), elec, vdw, net_force (`kcal/mol/Å`); meta lj_source,
charge_source, group_a, group_b, n_excluded_pairs.

**trace / proximity / atoms** (`analysis/structure.py`, internal). `structure_groups(project)`
(visual.structure_groups or all, known only), `structure_selection_text` (error `no residue groups
available for the structural view`), `structure_residues(ctx)` (atoms of the groups minus
`mg_selection` ions; sorted unique residues), `chain_ids` / `chain_names` (exist in core),
`residue_group_membership` / `residue_group_index` (failing group selections skipped),
`TracePoints` (first atom of `split_csv(trace_atoms)` or `["C4'","CA"]`, else heavy-atom centroid
— `ag.positions.mean(axis=0)` on **float32**), `LINK_MAX=9.0`, `guess_element` (structure.py's own,
not the selections one), `atom_bonds` (topology bonds among atoms, sorted unique pairs; else
first-frame distance guess, `pairs_within(xyz, xyz, 1.95, None)`, same/adjacent residue, no H-H,
H only <= 1.25).
- trace: settings visual.structure_groups, trace_atoms, align_frames, show_ions,
  mg_selection. Process: first frame = reference (+ links); later frames kabsch onto reference when
  align and >= 3 residues (ions transformed too). Arrays coords float32 [n, n_res, 3], ion_coords
  float32 [n, n_ions, 3], chain int64, group_index int64, group_member bool [n_groups, n_res],
  is_point bool, links int64 [m, 2], resindex, resid int64, resname `U8` (truncates at 8),
  chain_name `U*`. Meta residue_labels, ion_labels (`resname+resid`), group_names. Error
  `structural view: no residues selected`.
- proximity: `proximity_set` resolved (`resolve_selection` then `select("sel:...")`, heavy);
  member residues distance 0; struct heavy atoms of non-members; `pairs_within(struct, set, radius,
  box)` min per residue; `min_distance` float32 [n_res, n] initialised inf; arrays member bool,
  resindex; matrix min_distance kind "value" unit Å; meta residue_labels, radius.
- atoms: structure atoms (heavy unless rep_hydrogens), stable-sorted by index; coords float32
  [n, n_atoms, 3] superposed with the trace points' kabsch; row, bonds, resid, atom_index, resindex
  int64; meta atom_names, elements (`atom_elements`), resnames.

**basepairs** (internal). `engine` ("builtin"/"barnaba"), `stride`, `cutoff`, `wobble`,
`noncanonical`. Rows = `structure_residues` (error `base pairs: no residues selected for the
structural view`). Atoms per residue: heavy atoms, `*` -> `'`, first of a name wins. Builtin:
`secstruct.PairDetector(resnames, atom_index, cutoff, wobble).detect(xyz)` on frames `k % stride ==
0`, held between (`source[k]`). Barnaba: batch of 100 frames to a temp PDB -> `barnaba.annotate`.
Arrays map_* (value 1 WC, 2 wobble, 3 other), count int64, resindex int64; meta residue_labels,
engine, stride, value_labels, (barnaba) lw_first_frame.

### 1.3 Dependencies inside MDATK

| MDATK | Used by | MDBridge status |
|---|---|---|
| `spatial.pairs_within` (cKDTree open / `capped_distance` periodic), `minimum_image_vectors` (orthorhombic only, `np.round`, `np.allclose(angles, 90)`), `kabsch` (weighted, SVD, det correction, identity for < 3 points), `_usable_box` (shape 6, lengths > 0) | every geometric analysis | **missing** (periodic grid exists privately in `src/core/atom_selection.cpp`) |
| `secstruct`: `PAIR_RULES`, `C1_WINDOW`, `SCORE_TIE`, `PairDetector.candidates/detect`, `base_letter`, `canonical_name`, `lw_kind`, `dotbracket_to_pairs`, `pairs_to_dotbracket`, `_stems`, `unique_pairs`, `BarnabaRunner` | basepairs; panels, annotate (later) | **missing**; `parseResidueName` exists (`io/pdb.hpp`) |
| `data.py` AnalysisResult / Dataset / CSRBuilder | all | exists (`core/dataset.hpp`); NumArray lacks a string dtype (trace `resname`, `chain_name`) |
| `io/cache.py` CacheStore: `<cache>/<name>-<key>.npz` (`savez_compressed`, extra `_frames` int64, `_time_ps` float64) + `.json` (`project, mdatk_version, created, name, key, params, meta, arrays (sorted), n_frames`), manifest `MDATK_cache.json` (`project, mdatk_version, updated, system, frames{start,stop,step,n_frames,first_frame,last_frame}, time_axis{unit,n,start,end,dt=median(diff)}, groups, available=describe(), entries[{name,key,file,params,outputs sorted}], all_entries` merged with the old manifest, sorted), `write_json` (indent 2, `default`), `load_result`, `load_dataset_from_manifest` | pipeline | **missing** (own format allowed) |
| `pipeline.Pipeline`: `check` (validate + layout warnings), `resolve_inputs` (**cleaning runs cpptraj**), `universe`, `input_fingerprint` ({topology fp, trajectories fps, frames{start,stop,step}, dt_ps}), `required_analyses` (enabled non-internal + `render.panels.required_analyses` + chart series keys, minus "external", invalid instances skipped with `WARNING: skipping data source 'x': <exc>`, sorted by registry index then name), `analyze(names, force)` (per-name cache load + `"<label>: loaded from cache (<name>-<key>)"`, run missing in one pass, save, frames-only fallback, Dataset, `apply_display_overrides` (energy target), external data (`PipelineError("external data: ...")`), write manifest, progress `"analysis done"`), `load_cached_dataset` | `analyze` CLI | **missing**; validate, layoutWarnings, loadSystem, frameIndices, fileFingerprint, stableHash, loadExternal, Reporter exist |
| `render.panels.required_analyses` | required_analyses | per panel type: structure -> trace (+atoms if representations, +proximity if color_mode proximity or glow_mode != none); sequence -> trace (+proximity); secondary -> trace (+proximity) + basepairs unless `source == "dotbracket"`; others -> analysis part of `spec.series` keys |
| `cli.cmd_analyze` | `analyze -p ... [--only a,b] [--force]` + `_common` options; `pipe.check()`; prints `json.dumps({"cache_manifest", "n_frames", "data_keys": describe()}, indent=2)`; stderr via `make_reporter` / `ConsoleProgress` | **missing**; `_common` options exist for info/export-pdb |
| `cli.cmd_list` | prints analyses (`  %-10s label` + internal tag), conda envs with cpptraj (`cleaning.conda.find_cpptraj_envs`), panel types (`  %-12s label - description`), layout presets (`PRESET_LABELS`), fixed help text | **missing**; panel/preset tables and conda discovery absent in C++ |

---------------------------------------------------------------------------------------------------

## 2. External calls

### (a) MDAnalysis -> clean-room behaviour specs

| Call (as MDATK uses it) | Where | The spec must cover |
|---|---|---|
| `MDAnalysis.lib.distances.capped_distance(a.astype(float32), b.astype(float32), max_cutoff=c, box=dims (6,), return_distances=True)`; `min_cutoff=None`, `method=None` (auto) | `pairs_within` with a usable box (use_pbc only): contacts, hbonds (+H guess at 1.25), mg (2 calls), stacking (centroids vs centroids incl. self pairs i==i), energy, proximity | (1) pair set: inclusive/exclusive at the cutoff, in which precision the comparison is made (float32 vs float64), for each automatically chosen method and whether methods can disagree; (2) the method auto-selection rule (sizes, cutoff/box ratio) if outcomes differ; (3) returned distance values: dtype, precision (float32 arithmetic?), minimum image for orthorhombic and triclinic cells, points outside the primary cell / far outside, cutoff larger than half a box length (one pair per (i,j) or several images?), (4) order of returned pairs (sorted? method dependent) — matters for float sums and first-wins ties; (5) degenerate inputs reachable after MDATK's filter (`lengths > 0`): angles 0/180/NaN, flat cells, NaN coordinates, duplicated points, a and b the same array; (6) dtype/shape of `pairs` (intp [m,2]). Probe-vector generator included |
| `MDAnalysis.analysis.rms.rmsd(a, b, weights=None or masses float64, center=True, superposition=True)` on float64 (n,3) | rmsd without fit selection | centring (plain vs weighted mean with weights), weight normalisation, the optimal-rotation RMSD algorithm's numerics (QCP: eigenvalue iteration tolerance, max iterations, fallback when not converged), result for identical / near-identical structures (exactly 0 or ~1e-7 from cancellation?), 1 and 2 atoms, collinear/planar sets, zero total weight, NaN input, returned type; errors for mismatched shapes/weights (unreachable in MDATK; list only) |
| Universe / AtomGroup / Timestep behaviour used by the analyses beyond M2: `traj[int_array]` iteration (order, repeated indices, the current frame afterwards — rmsd seeks in setup), `ts.time`, `ts.dimensions` (None without a cell; zero cells), `hasattr(u.atoms, "charges")` per topology format, masses per format (parm7 / PDB guessed / GRO), `ag.bonds` membership (>= 1 atom in group) and `bonds.indices` order, `u.bonds` on bond-less formats (exception type), `atoms.positions` dtype, `ag.residues` / `residues.atoms` order, `atoms.fragindices`, `chainIDs` of first atoms, `ag[bool_mask]`, `ag | ag` | framework, all | Mostly already pinned by the M2 dumps (`tests/data/topology/*.topology.json`, `io/loader`, `io/trajectories`); the spec agent checks each item against those dumps and specifies only what is not covered (expected: fancy-index iteration + current frame, bonds-less exception, charges presence per format) |

`MDAnalysis.analysis.rms` and `capped_distance` are the only MDAnalysis *algorithms* the analyses
call; selections (`select_atoms`) are already reimplemented in `core/atom_selection.cpp`.

### (b) numpy / scipy (port with notice, or trivial)

| Call | Where | Decision / precision notes |
|---|---|---|
| `scipy.spatial.cKDTree(a64)` + `ta.sparse_distance_matrix(tb, r, output_type="ndarray")` | `pairs_within` without box (default `use_pbc=false`: the **main path**) | Pair set: `d <= r` (verify in scipy source); distance = sqrt of x,y,z squared sum in float64 from float32-widened coords. Output order = tree traversal order. Plan: own cell list (reuse/extract `CellGrid` from `atom_selection.cpp`) with a deterministic order, exact pair set; sums compared with tolerance. If tolerance fails, port scipy's `cKDTree` build (leafsize 16, sliding midpoint, balanced, compact) + `sparse_distance_matrix` traversal (BSD-3, with notice) to reproduce order exactly |
| `np.linalg.svd` (LAPACK gesdd) | `kabsch` (3x3), `ring_geometry` (batched k x 3, k = 5/6/9) | Reuse the one-sided Jacobi SVD in `src/core/rigid_transform.cpp` (expose in a small `core/linalg.hpp`). Sign ambiguity harmless (kabsch: U,V flip together; ring normal used via abs/squares). LAPACK differs in last bits -> tolerance (rigid cases already compare within a few 1e-15) |
| `np.linalg.det` | kabsch sign | trivial (3x3 cofactor); only the sign is used |
| `@` matmul (`w @ mobile`, `(x - mc) @ r.T + rc`, `((m-mc)*w).T @ (ref-rc)`) | kabsch, rmsd, trace, atoms | BLAS dot order differs -> 1e-16 relative; trace/atoms coords are then cast to **float32**: occasional 1-ulp differences; compare within 1 float32 ulp |
| `np.sum` / `np.mean` / `np.average` on 1-D float64 | rmsd (`sq.mean()`, average), energy `e_el.sum()`, `e_vdw.sum()` | numpy pairwise summation (port exists for float32 in `src/io/mrc.cpp`; generalise to double, same notice) — exact only if pair order matches |
| `sum(axis=0)` of (n,3) float64 | energy net force | sequential row accumulation (not pairwise) |
| `ag.positions.mean(axis=0)` on **float32** | TracePoints centroid | float32 accumulator, sequential over rows, float32 divide; reproduce exactly |
| `positions.mean(axis=1)` float64 | ring centroids | sequential over ring atoms |
| `np.linalg.norm(x, axis=-1)`, `(v1*v2).sum(1)` | hbonds, mg, stacking, PairDetector, links | `sqrt((x0²+x1²)+x2²)`, trivial |
| `np.round` (half-to-even) | `minimum_image_vectors`, PairDetector score grid `round(score/0.05)*0.05` | `std::nearbyint`; same IEEE ops in the same order |
| `np.cos(np.radians(x))` | hbonds, mg, stacking cut-offs | numpy may use SIMD/SVML kernels on AVX-512 (as found for `exp` in M2): the generator records the hex values numpy gives next to libm's; tests compare against libm |
| `np.unique(return_counts)`, `searchsorted`, `isin`, `add.at`, `minimum.at`, `repeat`, `cumsum`, `flatnonzero`, `triu_indices`, `split`, `argsort(kind="stable")`, `lexsort` (stable), `sort` | everywhere | trivial reimplementation; keep stable orders (mg ion grouping, PairDetector greedy order) |
| `np.median(np.diff(t))` | manifest dt | trivial |
| `np.savez_compressed` / `np.load` | cache | not reproduced byte-wise (cache incompatible by design); MDBridge writes its own files (proposal: `.npz` written with zlib so Python tests can read MDBridge's cache, or a plain binary) |
| `re.sub`, `re.findall` | `_out_key`, `_GROUP_TOKEN` | trivial; ASCII classes, runs -> one `_` (code-point vs byte irrelevant) |
| `stable_hash`, `file_fingerprint`, `json.dumps(indent=2)` | cache keys, CLI stdout | exist (`core/hash.hpp`, `core/json.hpp`) |

### (c) Already in MDBridge

`core::selectAtoms`/`select`/`resolveSelection`/`expandGroups`, `heavy`, `atomElements`,
`elementFilter`, `splitCsv`, `residueLabels`, `chainIds`/`chainNames`, `fragmentIndices`,
`AtomGroup` (`|`, masked, unique, residues, bonds), `Topology` (charges/masses/bonds optional),
`Frame` (float32 positions, optional box, time), `io::loadSystem`, `frameIndices`, `LjTable`,
`isParm7`, `loadExternal`, `core::Project` (`AnalysisParams`, `VisualParams`, `toPlain`, `validate`,
`layoutWarnings`, `cacheDir`, `groupNames`), `stableHash`, `fileFingerprint`, `Reporter`,
`Dataset`/`AnalysisResult`/`CsrBuilder`, `pythonJsonCompact`, `parseResidueName`, the rigid
transform SVD, numpy pairwise sum (float32), the CLI's `--set` overrides and common options.

---------------------------------------------------------------------------------------------------

## 3. `chem` functionality needed now

Only `basepairs` needs chem code. Create `src/chem` (module `mdbridge_chem`, depends on core, io)
with `secstruct.{hpp,cpp}` ported from `mdatk/secstruct.py`:

- **needed now:** `PAIR_RULES`, `C1_WINDOW`, `C1_IDEAL`, `SCORE_TIE`, `BASE_PAIR_ATOMS`,
  `base_letter`, `canonical_name` (both via `parseResidueName`), `PairDetector` (`candidates`
  in exactly MDATK's order: rule order, chunks of `2_000_000 // iq.size` rows, row-major
  `np.nonzero`; `detect`: stable lexsort by (rounded score, |C1'-C1' - 10.5|), greedy), `lw_kind`.
- **cheap, port now (pure, tested by `test_secstruct.py`, needed by the secondary panel and
  annotation next):** `dotbracket_to_pairs` (with its exact error texts), `_stems`, `_crosses`,
  `pairs_to_dotbracket`, `unique_pairs`, `KIND_LABELS`.
- **not ported:** `BarnabaRunner`/`import_barnaba` (GPL engine). `engine = "barnaba"` in MDBridge
  (owner decision): setup warns "basepairs: Barnaba is not available in MDBridge; using the builtin
  detector" and runs the builtin engine (params, and so the cache key, keep the project's values),
  recorded as an intended difference; parity cases skip barnaba. (MDBridge's own LW annotation is milestone 6.)
- structure helpers (`structure_residues`, TracePoints, `guess_element`, `atom_bonds`) stay in
  `analysis` (they are analysis code in MDATK).

---------------------------------------------------------------------------------------------------

## 4. Parity strategy

### 4.1 Systems (all generated, small, committed under `tests/data/analysis/`)

| System | Source | Exercises |
|---|---|---|
| `demo_dry` (parm7 + `io/trajectories/demo_dry.nc`, 24 frames; also `.xtc`, `.dcd`, `.trr` copies already in `io/trajectories/files/`) | existing | every analysis with MDATK's demo project (groups 1-12 / 13-24, energy eps 4r, target -60, Mg, H with bonds, charges, parm7 LJ); xtc = lossy coords, dcd/trr = other readers |
| `demo_dry` + variants of the project | generator (`--set`-style overrides) | use_pbc on (orthorhombic box of the demo), mass_weighted, fit_selection, superposition off, reference_frame (0, 5, 999 clamp), frames start/stop/step, contacts same group (A=A), heavy_only off, hbonds use_angle off, mg geometry on, bridge_shell outer, occupancy_selection, stacking scope inter, energy use_parm7 off (generic LJ), exclude_bonded 0/1/3, distance_dependent off, basepairs stride 1/10, wobble off, cutoff 3.0, trace align off / show_ions off / trace_atoms "P", proximity set `sel:resid 3`, atoms with hydrogens, instances `contacts@Target+Aptamer`, `rmsd@Aptamer`, errors (unknown group, group inside other) |
| `periodic_ortho`, `periodic_tric` | new generator: MDATK `write_parm7` + MDAnalysis NCDFWriter (like `make_trajectory_cases.py`); a small RNA fragment + Mg straddling cell faces, boxes 25 Å cubic and triclinic (60/70/80°), cutoff up to 12 Å vs box | the capped_distance path, minimum-image vectors (orthorhombic only by design), wrapped H guessing |
| `pdb_nobonds` | PDB of the demo frame 0 without CONECT, plus `.dcd`; `gro` variant | no charges (elec zero warning), generic LJ, H guessed from distance (bonds unknown), elements guessed, chainIDs |
| `noh` | demo heavy atoms only (PDB) | hbonds "no hydrogens" fallback |
| `aromatic` | hand-built PDB: PHE/TYR/TRP/HIS(HIE/HSD) next to bases, 3 frames | AROMATIC_RESIDUES rings, ring sizes 5/6/9 batches, planarity cut |
| MDATK micro systems (`tests/test_analyses.py`, `test_basepairs.py` helpers) | generator writes atoms + coordinates as JSON; C++ builds a `Topology` + `MemoryTrajectory` | each Python test's exact cases (contacts native fraction, H-bond angle linear/bent/far with and without bonds, Mg bridges/geometry/occupancy, stacking criteria, analytic Coulomb+LJ, energy bonded exclusions, rmsd superposition, base-pair rules / ties / stacked neighbours / L-RNA) |
| `demo_wet` (optional) | `build_demo(n_frames=6, solvent=True)` if < ~2 MB | Mg occupancy with water, larger neighbour searches, analyses without cleaning |

No real research data (README rule).

### 4.2 Expected data

`tests/parity/make_analysis_cases.py <MDATK checkout> tests/data/analysis` (run with
`/opt/mdbridge-venv/bin/python -I -B`, **`independent_tables.apply()` and
`residue_rule_universe`** before any Universe, since elements / masses / selections feed the
analyses): MDATK in-process, `load_universe` + `AnalysisContext` + `run_analyses` per (system,
project variant, analysis instance list). For each result: name, params (tagged JSON),
meta (JSON), `outputs`, every array as `{dtype, shape, values}` — integers as JSON, floats as
`float.hex()`, arrays > 4096 elements as SHA-256 of the little-endian bytes **plus** a sparse
sample (every 97th value) so tolerance comparison stays possible; times; Reporter log lines and
warnings; errors as `{type, message}`. Also: `cache_key` of each analysis for a *fixed* fingerprint
dict (keys independent of mtime), `expected_keys`, `key_label`, `referenced_groups`,
`instance_project` errors, `required_analyses` per project (panel sets). `pairs_within` /
`kabsch` / `minimum_image_vectors` vectors with random inputs (box and no box) + MDATK's
`secstruct` functions go to `analysis/spatial_cases.json` and `chem/secstruct_cases.json`.
Add each generator to `tests/parity/regenerate_all.py` (after the loader step) and document it in
`tests/data/README.md`. The spec agents' probe generators (`make_distance_cases.py`,
`make_rmsd_cases.py`) are separate steps (they call MDAnalysis; no MDAnalysis code is copied).

### 4.3 Tolerances

| Quantity | Comparison | Why |
|---|---|---|
| integer arrays (counts, CSR offsets/i/j, classes, links, bonds, indices, resids), bool arrays, labels, meta, outputs, dtypes, shapes, warnings, error texts, times | exact | discrete; times are read values |
| pair sets of `pairs_within` (no box) | exact (sorted) | same inclusive `<=` on the same float64 distances |
| pair sets with box | exact except pairs within ~1e-5 Å of the cutoff (flagged by the generator) | capped_distance precision unknown until the spec; selections already accept this (`atom_selection.hpp`) |
| any threshold decision (cutoffs, angle cos, planarity, C1' window, PairDetector score grid) | generator records the margin to the nearest threshold per frame/analysis; the C++ test accepts a different count only where the margin < 1e-9 (relative) | numpy/LAPACK last-bit differences cannot be ruled out |
| float64 analysis series (rmsd with own kabsch path / no superposition, q_native) | rel 1e-12 / abs 1e-12 | BLAS order; q_native is exact in practice |
| rmsd through `rms.rmsd` | abs 1e-6 Å near 0, rel 1e-9 otherwise (spec to confirm) | QCP iteration tolerance, sqrt of cancelled sums |
| energy total/elec/vdw | rel 1e-11 of the sum of |terms| (generator records it) | pair order of cKDTree vs ours, pairwise summation |
| net force | abs 1e-9 x sum of |force terms| | cancellation |
| float32 coords (trace, atoms, ion_coords), min_distance | exact for the reference frame and unaligned runs; aligned frames within 1 float32 ulp | float64 -> float32 rounding after BLAS-order differences |
| `cos(radians(x))` | exact vs libm; numpy's value recorded for information | SVML (see 2b) |

### 4.4 Cache and CLI parity

- Cache files: not compared (ARCHITECTURE). Tests check MDBridge reloads its own entries, that
  keys equal MDATK's for fixed inputs, that a parameter change invalidates exactly the affected
  analyses (MDATK `test_pipeline.py::test_full_run_and_cache_reuse`, `test_v002.py` keys), and that
  the manifest has MDATK's schema (`project`, `system`, `frames`, `time_axis`, `groups`,
  `available`, `entries`, `all_entries`; version field named MDBridge's).
- CLI: extend `tests/parity/cli_parity.py` with `analyze` (demo project, `--only` with instances,
  `--force`, `--frames`, `--set analyses.use_pbc=true`, unknown analysis -> MDATK traceback ->
  MDBridge `error: ...` exit 1, invalid project -> check() message, external data) comparing exit
  status, stdout JSON (normalised `{work}`), stderr lines (cache keys normalised `{key}` since
  fingerprints carry mtimes) and the set of files written (names, not bytes). `list`: everything
  except the conda section (environment dependent; compared with a fake conda root via env vars
  if the port supports it, else `exit_only` + the other sections line by line).
- Each Python test in `test_analyses.py`, `test_basepairs.py` (builtin cases), `test_secstruct.py`,
  the analysis parts of `test_pipeline.py` / `test_v002.py` gets a C++ counterpart
  (`tests/analysis/test_*.cpp`, `tests/chem/test_secstruct.cpp`).

---------------------------------------------------------------------------------------------------

## 5. Work breakdown

Notation: **[S]** clean-room spec (reads MDAnalysis, writes only a behaviour description + a probe
generator that *calls* MDAnalysis; never writes MDBridge code), **[I]** clean-room implementation
(reads only the spec + MDBridge; never MDAnalysis source/docs beyond the spec), **[P]** MDATK port
(may read MDATK, numpy/scipy; never MDAnalysis), **[T]** test-data generator, **[C]** CLI. Critical
path marked ★. Specs live in `docs/specs/m3-*.md` (decision for the owner: commit them, or keep
them outside the repo as for M2).

### Wave 0 (all independent, start together)

**W0.1 [P] ★ Analysis framework** — `src/analysis/` (`CMakeLists.txt`, `include/mdbridge/analysis/
analysis.hpp`, `registry.cpp`, `context.cpp`, `runner.cpp`, `common.cpp`, `keys.cpp`),
`tests/analysis/test_framework.cpp`. Contents: `Analysis` interface (`setup(frames)`,
`process(k, const Frame&)`, `finalize()`), static descriptors (name, label, version, internal,
groupSlots, outputLabels, enabled, paramsFor -> `core::Value`, selectionTexts, expectedOutputs,
withGroups), registry in MDATK order, instance parsing/errors, `instanceProject` (copy),
`keyLabel`, `expectedKeys`, `referencedGroups`, `cacheKey`, `AnalysisContext`, `runAnalyses`
(single pass; frame-parallel later), `residueRows`, `without`, `perRowCounts`; `pairsWithin`
**interface** with a brute-force reference implementation (lets analysis ports start before W1.1).
Also a small core change: `NumArray` string dtype (`<U n>`, for trace) — coordinate with dataset
tests. Deps: none. Accept: framework cases of W0.6 (keys, labels, errors) exact.

**W0.2 [S] ★ capped_distance spec** — `docs/specs/m3-capped-distance.md`,
`tests/parity/make_distance_cases.py` -> `tests/data/analysis/distance_cases.json` (random and
crafted inputs: boundary distances in float32, ortho/triclinic, outside-cell points, big cutoffs,
self arrays, every auto method region). Covers the list in 2(a). Deps: none.

**W0.3 [S] rms.rmsd spec** — `docs/specs/m3-rmsd.md`, `tests/parity/make_rmsd_cases.py` ->
`tests/data/analysis/rmsd_cases.json`. Covers 2(a). Deps: none.

**W0.4 [S] Universe-behaviour gap check** — short `docs/specs/m3-universe-notes.md`: which items of
2(a) row 3 are already pinned by M2 data, specs + probe cases for the rest. Deps: none.

**W0.5 [P] chem/secstruct** — `src/chem/` (new module: CMake, `include/mdbridge/chem/secstruct.hpp`,
`secstruct.cpp`), `tests/chem/test_secstruct.cpp`, generator `tests/parity/make_secstruct_cases.py`
-> `tests/data/chem/secstruct_cases.json` (dot-bracket parse/write/errors incl. 26 letter levels,
`unique_pairs`, PairDetector candidates and detect on `test_basepairs.py` geometries, perturbed and
random frames, score ties at the 0.05 grid, missing atoms, L-RNA names, every residue-name scheme).
Deps: none (uses `parseResidueName`). Accept: exact (float columns of candidates rel 1e-15).

**W0.6 [T] ★ analysis test-data generator** — `tests/parity/make_analysis_cases.py` (+ system
builders: periodic, pdb_nobonds, noh, aromatic, micro systems), `tests/data/analysis/` per 4.1/4.2,
`regenerate_all.py` step, README section. Deps: none (MDATK + MDAnalysis venv). Deliver in two
drops: (a) framework cases + micro systems + demo_dry default project (unblocks W0.1/W2.*),
(b) all variants, periodic systems, margins. Accept: `regenerate_all.py --only analysis`
reproducible byte for byte.

**W0.7 [P] spatial (non-PBC) + kabsch** — `src/analysis/spatial.{hpp,cpp}`: `pairsWithin` open
path (cell list; extract `CellGrid`/`Geometry` from `atom_selection.cpp` into
`core/detail/neighbors.hpp` if suitable), `minimumImageVectors` (MDATK semantics: orthorhombic only,
`allclose(angles, 90)` with numpy's rtol 1e-5 / atol 1e-8), `usableBox`, `kabsch` (uses the
Jacobi SVD exposed as `core/linalg.hpp`), float64 pairwise sum (generalise the numpy port in
`mrc.cpp`; notice). Tests `tests/analysis/test_spatial.cpp` vs `spatial_cases.json` and MDATK's
`test_pairs_within_*` / `test_kabsch_*`. Deps: W0.1 interface. Accept: exact pair sets, distances
exact, kabsch rel 1e-14.

### Wave 1 (consume specs)

**W1.1 [I] ★ periodic neighbour search** — implements `pairsWithin` with a box from
`m3-capped-distance.md` only (may reuse MDBridge's own `CellGrid`, which is not MDAnalysis code).
Tests vs `distance_cases.json`. Deps: W0.2, W0.7. Accept: per the spec's stated tolerance; pair
sets exact away from the boundary band.

**W1.2 [I] optimal-superposition RMSD** — `src/analysis/rmsd_qcp.{hpp,cpp}` from `m3-rmsd.md`
(and the published QCP papers / Theobald's BSD reference implementation if the spec allows; not
MDAnalysis' translation). Tests vs `rmsd_cases.json`. Deps: W0.3. Accept: spec tolerance.

### Wave 2 (MDATK ports of the analyses; each starts on the brute-force `pairsWithin`, final
acceptance after W1.1 for the PBC variants)

| Chunk | Files (`src/analysis/`, `tests/analysis/`) | Deps | Acceptance |
|---|---|---|---|
| **W2.1 [P] contacts + hbonds** | `contacts.cpp`, `hbonds.cpp`, tests | W0.1, W0.7, W0.6a | MDATK micro tests + demo variants exact (counts, CSR, q_native, per_res) |
| **W2.2 [P] mg + stacking** | `mg.cpp`, `stacking.cpp` (batched ring SVD via `core/linalg`), tests | W0.1, W0.7, W0.6a | exact counts / classes; margins rule for planarity / angles |
| **W2.3 [P] energy** | `energy.cpp` (GENERIC_LJ, bonded exclusions BFS, LjTable), tests | W0.1, W0.7, W0.6a | analytic test exact-ish (rel 1e-14), demo rel 1e-11; exclusions exact |
| **W2.4 [P] rmsd** | `rmsd.cpp` (`_out_key`, mass weighting, fit selection, NaN rows, reference seek) | W0.1, W0.7, W1.2 | per 4.3 |
| **W2.5 [P] ★ structure: trace, proximity, atoms** | `structure.cpp` (structure groups/residues, membership, TracePoints with float32 centroid, guess_element, atom_bonds, links) | W0.1, W0.7, W0.6a | exact topology arrays; float32 coords per 4.3 |
| **W2.6 [P] basepairs** | `basepairs.cpp` (builtin engine, stride hold, barnaba -> builtin engine with a warning (owner decision)) | W0.5, W2.5 (`structureResidues`) | exact on demo (stride 1 and 10) and micro cases |

### Wave 3 (pipeline and CLI)

**W3.1 [P] ★ cache + analyze pipeline** — `src/analysis/cache.{hpp,cpp}` (store, manifest,
reload; own file format), `src/analysis/pipeline.{hpp,cpp}` (`requiredAnalyses` incl. the panel
table of 1.3, `inputFingerprint`, `analyze(names, force, reporter)`, display overrides, external
data, `loadCachedDataset`), cleaning-enabled projects: decision needed (see risks) — proposal:
use the cleaned files when MDATK's cleaning output for the same key exists, else error
`cleaning is not available in MDBridge yet`. Tests: `test_pipeline.cpp` (cache reuse and
invalidation, frame selection + external, invalid parameters, cancellation). Deps: all W2 for full
acceptance (can start on W0.1 + any one analysis).

**W3.2 [C] ★ `analyze` and `list` commands** — `src/app/commands/commands.cpp` (+ `analyze.cpp`,
`list.cpp`), `ConsoleProgress` / reporter output as MDATK (`[ 12.3%] msg`, tty handling),
`tests/cli/test_commands.cpp` cases, `cli_parity.py` cases (4.4). `list` needs static tables:
panel types (name, label, description) and layout presets — generate with a small
`tools/gen_panel_tables.py` from MDATK (data only, MDATK's own) into `src/core/panel_tables.inc`
(render milestone reuses them); conda env discovery (`find_cpptraj_envs`) ported from
`mdatk/cleaning/conda.py` or stubbed (decision). Deps: W3.1.

### Wave 4 (after parity is green)

**W4.1 [P] frame-parallel runner** — split `process` into a sequential reader + per-frame workers
(analyses keep per-frame outputs indexed by k; stateful first-frame steps — rmsd reference,
hbonds H guess, trace/atoms reference + links, basepairs stride — run in a `prepare(first frame)`
phase); results identical for any thread count (test runs 1 and N threads, byte-equal arrays).
Deps: W2.*, W3.1.

**W4.2 docs** — ARCHITECTURE status, LICENSING / THIRD_PARTY_NOTICES rows (numpy pairwise double,
any scipy port), README of test data.

### Order and critical path

```
W0.2 spec ──► W1.1 PBC search ──────────────┐
W0.3 spec ──► W1.2 QCP ──► W2.4 rmsd ───────┤
W0.1 framework ─► W0.7 spatial ─► W2.1..W2.3│
                               └► W2.5 ─► W2.6 (needs W0.5)
W0.6 generator (a ► b) ─────────────────────┴─► W3.1 pipeline ─► W3.2 CLI ─► W4.1
```
★ critical path: **W0.1 → W0.7 → W2.5 → W2.6 → W3.1 → W3.2**, gated in parallel by
**W0.2 → W1.1** (all `use_pbc` parity cases) and **W0.6** (every acceptance test). Maximum useful
parallelism: 7 agents in wave 0, 6 in wave 2.

---------------------------------------------------------------------------------------------------

## 6. Risks and open decisions

1. **capped_distance numerics** (precision of distances/comparison, boundary, pair order, method
   switching) are unknown until W0.2; if they are float32 internally, periodic results will need
   the boundary-band tolerance permanently (as selections already do).
2. **Pair order → float sums** (energy, force) and **first-wins ties** (hbonds H guessing on
   symmetric geometry): exact parity would require porting scipy's cKDTree traversal (allowed,
   BSD) and matching capped_distance's order from the spec; plan accepts tolerance first.
3. **Threshold flips** from last-bit differences (LAPACK SVD planarity, BLAS kabsch, numpy SVML
   `cos`, PairDetector's `round(score/0.05)` grid on the noisy demo) — handled by the generator's
   margin records; a flip on the demo would surface as a count difference to be explained.
4. **QCP vs SVD**: rmsd of the reference frame (≈0) and nearly identical frames depend on the QCP
   iteration; tolerance must come from the spec, not guessed.
5. **float32 details**: TracePoints centroid in float32, trace/atoms coords cast to float32 after
   double math; 1-ulp differences expected on aligned frames.
6. **Barnaba engine**: unavailable (GPL) — intended difference; projects with `engine="barnaba"`
   fall back to the builtin engine with a warning (owner decision).
7. **Cleaning in `analyze`**: MDATK runs cpptraj when `cleaning.enabled`; MDBridge has no cleaning
   yet. Decide: error, reuse an existing cleaned output, or run cpptraj (scripts are a later
   milestone).
8. **`list`** depends on render tables (panel types, presets) and conda discovery, neither in C++
   yet; small generated table + port/stub needed.
9. **NumArray** has no string dtype (trace `resname` U8 truncation, `chain_name`) — core change.
10. **Cache format** is MDBridge's own; keys kept MDATK-compatible. Confirm with the owner that
    reading MDATK caches is not wanted.
11. **Clean-room hygiene**: the existing periodic grid in `atom_selection.cpp` may be reused by W1.1
    only if it was itself written clean-room (it states "computed independently of MDAnalysis");
    spec agents must not touch `src/`, implementers must not open MDAnalysis.
12. **Performance**: brute-force `PairDetector` candidates are O(n_p n_q) per rule like MDATK;
    fine for test systems, revisit with W4.1.
