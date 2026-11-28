# Milestone 7 plan: the `workflow` module

Goal (ARCHITECTURE "Migration order" 7 + "Definition of done"): MDATK v0.0.4.1's workflow layer in
MDBridge: the sixteen step types of the Project Planner with their checks and script generators, the
planner (plan model, resolution, file naming, master script and summary, templates), the stand-alone
script wizard, the runner that starts the generated scripts on request, the `workflow` command group
(`tools`, `check`, `write`, `run`, `generate`, `template`, `restraints`) and the GUI parts (planner
canvas, wizard). The generators are written in Lua and loaded from `.lua` files, so users can add step
types (ARCHITECTURE "Modules", `workflow` row; "Scripting: embedded Lua for script/wizard
generators"). The GUI parts become the second program `rnado` (owner decision D1 of milestone 5).

Licensing: MDATK code (owner's) is ported freely. The workflow layer does not use MDAnalysis
(section 3.1), so no clean-room specs are needed. New third-party material is the Lua runtime and a
few published parameter tables (section 3); the external programs the scripts drive are behaviour to
document, never code to copy.

References: MDATK v0.0.4.1 = `htchamblee/MDATK` tag `v0.0.4.1` (9875bb6, the parity reference);
"0.1.2" = its later `main` (2d93338, version 0.1.2), read for sections 2 and 5. Line numbers are
v0.0.4.1's unless a version is named. MDBridge base: a977c7c (milestone 6 in progress: W2-A and W3-C
open, `docs/specs/m6-status.md`). "Ran" marks behaviour confirmed by running MDATK 0.0.4.1 (Python
3.13, the scratch venv; section 1.11), "read" behaviour taken from the code.

---------------------------------------------------------------------------------------------------

## 1. Inventory of MDATK 0.0.4.1's workflow layer

Sizes (lines): `workflows/` 8 684 in 25 files: `amber.py` 1 607, `refine.py` 1 000, `plan.py` 618,
`farfar2.py` 480, `mg_params.py` 435, `gromacs.py` 429, `dictstep.py` 362, `base.py` 335,
`buildstep.py` 322, `simrna.py` 310, `vmd.py` 293, `rna.py` 272, `schema.py` 269, `rnamake.py` 205,
`runner.py` 187, `misc.py` 97, `naming.py` 73, `registry.py` 44; the embedded helpers `scripts/` 843
(section 1.4). GUI: `gui/planner.py` 1 056, `gui/wizard.py` 910 (GTK 3), hooks in `gui/app.py`.
Command line: `cli.py:371-417, 574-760, 934-994`. None of it imports MDAnalysis, NumPy or SciPy
(grep, read); the generators call the chem layer (m6-plan 4.4) and `citations`.

"MDBridge" columns: **done** = exists with parity tests (milestone and id), **part** = some of it,
**-** = missing.

### 1.1 What MDBridge already has

| Piece | MDATK | MDBridge |
|---|---|---|
| Step-type schemas: `Option` (8 kinds, `coerce`, `requires`), `Section` (repeated records, templates), `Port`, `Issue`, `FILE_KINDS` (17), `TOOLS`, `CATEGORIES`, `normalise`, `defaults` | schema.py:21-269, registry.py:11-44, base.py:144-225 | done (M5: `core::StepType`, src/core/include/mdbridge/core/step_types.hpp, generated from MDATK by tools/gen_step_types.py into src/core/step_types.inc; tests/data/step-types) |
| File naming `<project>_<step>_<tool>_<stamp>.<ext>`, step folders, `parse_name` | naming.py:26-73 | done (M5: core/naming.hpp) |
| Plan model `PlanNode`, `PlanEdge`, `WorkflowParams` (+ 0.1.x keys kept), path resolution of `files` and `path_settings`, comments of every setting | config.py:414-461, 540-543; workflows/__init__.py:24-33; projectfile.py:200-260 | done (M5: core/project.hpp:642-733, PRJ-09; tests/data/project-file/workflow) |
| Step log, `logged`, output tail (20 000 characters) | projectfile.py:51, 389-480 | done (M5: core/step_log.hpp) |
| Citation block of a script header (`comment_lines`, `step_resources`) | citations.py:61-142 | done (M5: core/citations.hpp:66-88) |
| Process runner (PATH lookup, merged output, cancel ending the process tree) | runner.py:56-161 (subprocess) | part (M5: core/process.hpp, used by `clean`; CLN-03) |
| Trajectory cleaning (cpptraj script, `conda run`, reuse) | cleaning/cpptraj.py, conda.py | done (M5 W0-E: `clean`; DIFFERENCES CMD-01 ... CMD-11) |
| The chem functions the generators call | m6-plan 4.4 | done (M6: chem/generator_api.hpp) |
| Workflow menu entries (planner, wizard: disabled shells; PDB clean-up, builder: windows) | app.py:158-168, 1595-1610 | part (M5 D2, M6 W3-B: src/app/gui/main_window.cpp:386-453, gui/tools/host.hpp) |
| Plan resolution, checks, generators, writer, runner, `workflow` commands, planner, wizard | | - (commands.hpp: "MDATK's workflow commands follow with milestone 7") |

### 1.2 The generator framework (`base.py`, `schema.py`, `plan.py`)

| Part | MDATK (file:line) | What it does |
|---|---|---|
| `PortValue` | base.py:21-39 | files along a connection: `paths` (plan-relative, or absolute when hand-picked), `kind`, `data` (facts known at plan time: `rna_model`, `chirality`, `nbfix`, `hmr`, `c4`, `mg_pairs`, `velocities`, `box`, `pairs` ...), `source`, `source_label` |
| `StepContext` | base.py:42-120 | names (`name`, `own`, `from_own`), inputs (`has`, `value`, `inp`, `inps`, `data`, `rel` = relative to the step folder by `posixpath.relpath`, `abspath`, `existing` = the absolute path when the file exists now, `source_label`), `groups`, `preview` (the wizard's preview: no network or slow work) |
| `Tool` | base.py:123-243 | `name`, `label`, `category`, `tag`, `description`, `generates`, `sections`, `inputs` / `outputs` (or `input_ports(s)` / `output_ports(s)` when settings change them), `path_settings`; `check`, `output_values`, `generate`, `describe` (box text), `resources` (citation ids) |
| `normalise` | base.py:168-211 | defaults overlaid with coerced settings; bad values -> default + warning issue "<message> - using the default"; unknown keys logged and dropped (MDBridge keeps them, PRJ-09) |
| Shell helpers | base.py:251-335 | `env_section` / `env_lines` (conda activation, setup lines), `script_header` (shebang, title, "Project P, step NN - written by MDATK 0.0.4.1 on <date>", "MDATK only writes this script; review it, then run it with: bash <script>", `set -eo pipefail`, `cd "$(dirname "$0")"`), `require_file` (`[ -s f ] || { echo "ERROR: ..." >&2; exit 1; }`), `fortran_float` (10.0, 0.002, 1.0e-6), `shell_quote`, `extra_lines`, `extra_flags` (`#` comments, also after a flag) |
| `embed` | scripts/__init__.py:15-37 | a helper module's source (`inspect.getsource`) without its docstring and `__future__` line, a new docstring, `# --- settings written by MDATK ---` with constants as Python `repr`, and `sys.exit(<call>)` |
| `Plan` resolve | plan.py:150-318 | duplicate ids, unknown types, invalid edges (warnings), topological order keyed by (step, y, x, id) (:212-234: independent steps follow the canvas position), per step: inputs from connections (merged for `multiple` ports) or hand-picked files (`;`-separated, `abspath(expanduser)`, :268-278), the "connect a step that produces it or choose a file" error, `check` and `output_values` (an exception becomes "internal error while checking this step: <text>", :286-291), (step, tag) collisions and shared step numbers (:296-317) |
| `Plan.write` | plan.py:348-382 | steps with errors skipped, a file written by two steps skips the second (:359-365), atomic writes `<path>.tmp` + `os.replace`, mode 0755 for `run`, `run-separate`, `helper` (:471-478), the stamp fixed on first write (:374-375), master `<project>_00_plan_<stamp>.sh` and summary `.md` only for a whole-plan write |
| Master script, summary | plan.py:384-448 | master: `bash "<rel>"` per `run` file in plan order, `#   run separately: bash <rel>`, `#   (not written: the step has errors)`; summary: folder, type, inputs (with source), outputs, files written, checks (no info lines), notes as quotes |
| `generate_step` | plan.py:451-468 | inserts the "Please cite" block after the header comments of `run` / `run-separate` scripts |
| Stand-alone wizard, analysis link | plan.py:486-517 | `scratch_node` (one settings record per step type in `workflow.scratch`), `single_step_plan`, `analysis_inputs` |
| Templates | plan.py:525-618 | `amber`, `farfar2-amber`, `builder-amber`, `cryoem-refine`, `amber-gromacs`; columns 280 px, rows 190 px |
| Graph edits | plan.py:71-147 | `new_node_id` (`n<k>`), `next_step`, `add_node`, `remove_node`, `connect` (kind check, no loops, one edge into a single-file input) |

### 1.3 Step types (16; 357 options)

Inputs -> outputs are the ports with default settings (`workflow tools`, ran); options counted over
all sections (ran). Files are those written with default settings and the variants of section 1.11
(`NN` = a two-digit number, `<stage>` = one per MD stage); every name is `<project>_<step>_<tag>_<stamp>.<ext>`.

| Step type (file:line) | Category | Ports | Options | Files written (role) | Programs the scripts run | Plan-time reads, chem calls |
|---|---|---|---|---|---|---|
| `amber-build` (amber.py:252) | Amber | structure (pdb) -> topology, coordinates, pdb (+ dry_topology, dry_pdb with `two_pass`; :406) | 36 | `amber-build.sh` (run); `tleap.leap`, `parmed.in`, `tleap-dry.leap`, `parmed-dry.in`, `mg.frcmod` (input); `prep.py`, `pdbtools.py`, `cp-lib.py`, `tip4p2005.py` (helper) | tleap, ParmEd, ambpdb, python3 | input PDB: L-RNA detection (`parse_ranges`, `read_pdb`, `assign_chirality`, :533-577), solute residue count (:228-245) |
| `amber-restraints` (amber.py:966) | Amber | topology (parm7), secstruct (opt.) -> restraints (dist) | 12 | `amber-restraints.sh` (run), `helix.py` (helper), `helix.dist` (input, when the topology exists) | python3 | a secondary-structure file (:155-201); the parm7 when it exists (:992, 1031-1033) |
| `amber-md` (amber.py:1186) | Amber | topology, coordinates, restraints (opt.) -> topology, restart, trajectory (multiple), dry_pdb (+ dry_topology, dry_trajectory; :1253) | 57 (10 + 41 per stage + 4 + 2) | `amber-md.sh` (run), `amber-prod.sh` (run-separate), `<stage>.in`, `<stage>-cont.in` (input) | pmemd, pmemd.cuda, pmemd.MPI (`mpirun -np {n}`), sander, cpptraj | the coordinates' PDB for the solute mask |
| `structure-builder` (buildstep.py:30) | Structure building | pdb, map (by mode) -> pdb, secstruct | 26 | `build.sh` (run), `helix.secstruct` (input) | `python3 -m mdatk build helix / mutate / fit`, or the AppImage (`$APPIMAGE`, :259-268) | `resolve_partner`, `resolve_third`, `pair_types`, `check_helix`, `parse_mutations`, `read_pdb`, `helix_kind` |
| `restraint-dicts` (dictstep.py:92) | Refinement | structure (opt.) -> dictionaries (cif, multiple) | 14 | `dicts.sh` (run), `dictionaries.py` (helper: MDATK's whole `dictionaries.py`, :330-335), `<code>.smi` (input) | acedrg or phenix.elbow, python3 | model codes (`read_pdb`, `assign_chirality`, `convert_names(ccd)`), `COMPONENTS`, `check_smiles`, `parse_smiles`, PubChem only when asked and not previewing |
| `farfar2` (farfar2.py:30) | Structure prediction | structure, map (opt.) -> silent (multiple), fasta, secstruct | 38 | `farfar2.sh` (run), `farfar2.flags`, `farfar2.fasta`, `farfar2.secstruct` (input) | rna_denovo (xargs -P) | `rna.build_model` |
| `farfar2-models` (farfar2.py:362) | Structure prediction | silent (multiple) -> best, models (multiple) | 14 | `farfar2-models.sh` (run) | rna_cluster, extract_pdbs, x3dna-dssr | - |
| `gromacs` (gromacs.py:88) | GROMACS | topology (parm7 / top), coordinates (rst7 / gro) -> topology, structure, trajectory | 37 | `gromacs.sh` (run), `gromacs-prod.sh` (run-separate), `<stage>.mdp` (input), `convert.py` (ParmEd Amber -> GROMACS, keeps NBfix / Mg pair terms), `posre.py` (helper) | gmx grompp / mdrun, ParmEd or acpype | port data (`nbfix`, `c4`, `hmr`, `mg_pairs`) |
| `phenix-rsr` (refine.py:597) | Refinement | model, map, dictionaries (opt.), secstruct (opt.) -> refined, log | 23 | `phenix-rsr.sh` (run), `basepairs.eff` (input) | phenix.real_space_refine (1.x / 2.x) | REMARK 950 PAIR / STACK (`read_annotation`), `read_pdb`, `_nucleotide`, `assign_chirality`, `lw_kind` |
| `coot-rsr` (refine.py:863) | Refinement | model, map, dictionaries (opt.) -> refined, log | 11 | `coot-rsr.sh` (run), `coot-rsr.py` (input: MDATK's Coot script, 0.9 and 1.x) | coot --no-graphics | residue names of the model |
| `rnamake` (rnamake.py:20) | Structure prediction | structure -> designs (multiple), best | 27 | `rnamake.sh` (run) | design_rna_scaffold | - |
| `simrna` (simrna.py:39) | Structure prediction | structure, secstruct (opt.) -> best, models, trajectory | 31 | `simrna.sh`, `simrna-cluster.sh` (run), `simrna.config`, `.seq`, `.ss`, restraints (input) | SimRNA (xargs -P), clustering, SimRNA_trafl2pdbs | `rna.build_model` |
| `vmd` (vmd.py:66) | Visualisation | topology, trajectory (opt.) -> image | 27 | `vmd.sh` (run), `vmd.tcl` (input) | vmd, ffmpeg (movie) | residue groups (`selections.parse_residue_ranges`) |
| `note`, `input-file`, `mdatk` (misc.py:15, 31, 69) | Planning, Analysis | -; - -> file (kind chosen); topology, trajectory -> - | 1, 2, 1 | none (`generates = False`) | - | `input-file`: existence check |

Plan-level files: `<project>_00_plan_<stamp>.sh` (master) and `.md` (summary). File kinds in all:
17 bash scripts (16 per-step, the master), 8 embedded Python helpers, about 15 program-input formats
(tleap, ParmEd, frcmod, pmemd mdin, DISANG, SMILES, rna_denovo flags, FASTA, dot-bracket, GROMACS mdp,
Phenix PHIL, Coot Python, SimRNA config / sequence / secondary structure / restraints, VMD Tcl) and the
Markdown summary.

### 1.4 Embedded helper scripts (`embed`, standard library only)

| Helper | Source (lines) | Written by | What it does |
|---|---|---|---|
| `helix.py` | scripts/helix_dist.py (102) | amber-restraints | DISANG restraints from a parm7 (atom indices by residue number and name); constants `PAIRS`, `BONDS`, `WINDOW`, `OFFSET` |
| `prep.py` | scripts/prep_pdb.py (104) | amber-build | hydrogens off, 5'-phosphates off, TER at numbering gaps, waters off, first model |
| `cp-lib.py` | scripts/cp_retype.py (160) | amber-build (`cp_source = lib`) | RNA_CP.lib and frcmod.cp from AmberTools' RNA.lib / parm10.dat at run time; CP values of Steinbrecher et al. 2012 (:27-29) |
| `tip4p2005.py` | scripts/tip4p2005.py (182) | amber-build (`water = tip4p2005`) | frcmod.tip4p2005 and tip4p2005box.off from AmberTools' TIP4P-Ew files (Abascal & Vega 2005 constants) |
| `convert.py` | scripts/amber2gmx.py (172) | gromacs | ParmEd conversion that keeps off-diagonal LJ pairs (needs ParmEd on the run machine) |
| `posre.py` | scripts/gmx_posre.py (86) | gromacs | POSRES includes per molecule |
| `pdbtools.py` | MDATK's `pdbtools.py` (1 210 lines written) | amber-build (`lrna`, default on) | `clean --scheme amber [--l-ranges]` before tleap (amber.py:849-858) |
| `dictionaries.py` | MDATK's `dictionaries.py` (1 448 lines written) | restraint-dicts | `finish` (renaming to CCD atom names) after acedrg / eLBOW |

The last two carry MDATK 0.0.4.1's chem code, including the bugs milestone 6 fixed (CH-01 ... CH-15)
and MDATK's dictionary table that disagrees with the CCD (CHD-01, guanine H21 / H22). Section 4.6 and
D3 decide what MDBridge writes instead.

### 1.5 Planner and project model

The plan is `[workflow]` in the project file (config.py:445-461): `project_name`, `stamp` (fixed on
first write), `use_stamp`, `out_dir` (default `<output.out_dir>/<project_name>_workflow`,
config.py:559-565), `step_folders`, `nodes`, `edges`, `scratch`. A node: `id`, `tool`, `step`, `title`,
`x`, `y`, `settings`, `files` (config.py:414-431). MDBridge reads and writes all of this already
(1.1). Only `use_stamp`, `step_folders` and the step numbers change file names; titles and positions
change the master order only through the topological key (plan.py:219-221, read).

### 1.6 Runner (`runner.py`)

`plan_jobs` collects the `run` files of the steps without errors in plan order (:88-107);
`run_jobs` runs them one after the other with `bash <path>` in the script's folder, a new session
(process group), `frozen.child_env()` (frozen.py:108), output streamed line by line, one step-log entry
per script (`step` "NN", the step's title and citations, `command`, `folder`, `output`, `summary`
"exit status N"), the first failure stops (:164-187). Cancel: SIGTERM to the group, SIGKILL after 10 s
(:56-69). Windows: never runs, with the message of `runner_problem` (:42-53). `separate_scripts` lists
the `run-separate` files, which are never started (:110-118).

### 1.7 Cleaning

`cleaning/cpptraj.py` and `conda.py` (the trajectory cleaning, not a planner step) are done in
milestone 5 (1.1). M7 reuses the process runner and `conda run` handling for the planner's scripts.

### 1.8 Command line

All options and outputs read from cli.py:934-994 and confirmed by running (section 1.11). Errors print
`error: <message>`, exit 1 (cli.py:1014-1019); usage errors exit 2.

| Command (cli.py) | Options | Output, step log | MDBridge |
|---|---|---|---|
| `workflow tools [TOOL]` (:574-600) | | categories and step types, or one step type's ports and options (`key kind {choices} default=<repr>` + tooltip) | - |
| `workflow check` (:616-623) | naming options: `--project/-p/--params`, `--out-dir/-o`, `--project-name`, `--stamp`, `--no-stamp`, `--quiet/-q` | issues on stderr (`LEVEL    label: message`), JSON `{steps, errors, folder}`; exit 1 with errors | - |
| `workflow write` (:626-654) | naming + `--step N`, `--force`, `--no-save`, `--run` | JSON `{folder, written, skipped, master, summary}`; step log `workflow write`, "Generated the plan's files"; the new stamp saved into the project (`_edit_project_workflow`, :737-743) | - |
| `workflow run` (:657-684) | naming + `--step N` | script output; one entry per script; `not started (run separately): <path>` | - |
| `workflow generate TOOL` (:686-715) | naming + `--step` (1), `--set/-s KEY=VALUE` (records `section.index.key`, :375-406), `--input/-i PORT=PATH`, `--dry-run`, `--run` | `===== <name>` + content (dry run) or JSON `{folder, written, stamp}` | - |
| `workflow template [NAME]` (:718-734) | `--project` (saved in place) | the list, or JSON `{added, project}` | - |
| `workflow restraints` (:746-760) | `--top`, `--sequence`, `--dotbracket`, `--offset`, `--wobble`, `--r R1 R2 R3 R4`, `--rk RK2 RK3`, `--output/-o` (helix.dist) | the DISANG file; issues on stderr; no step log | - |

### 1.9 GUI (GTK 3 in 0.0.4.1)

| Element | MDATK | MDBridge |
|---|---|---|
| Project Planner window: canvas (boxes with ports, connections drawn as curves, drag, rubber-band connect, Delete, Ctrl+wheel zoom, context menu: wizard, duplicate, delete, use for analysis), toolbar (*Add step* by category, *Templates*, *Delete*, *Arrange*, zoom, *Fit*, *Save project*, *Write all files*, *Run…*, *Open folder*), naming bar (project, time stamp, *New stamp*, step folders, folder), details pane (title, step number, issues, ports, files), status line | gui/planner.py:67-1056 (toolbar :610-664, naming :666-700, details :702-830, write / run :971-1003, use for analysis :1012-1022) | - (disabled shell, Ctrl+P reserved) |
| Script wizard: tool list, forms built from the schema (`OptionForm`: check boxes, radio buttons, combo boxes, spin buttons, text, multi-line, file fields; greyed by `requires`), repeated sections (`RecordEditor`: list, *Add* from templates, duplicate, remove, move), files page (connected / chosen inputs), issues list (click selects the field), live preview of every file, *Write*, *Write and run*, *Defaults*, *Open folder*; stand-alone (scratch settings) or a box's *Edit in wizard…* | gui/wizard.py:44-910 | - (disabled shell, Ctrl+Shift+W reserved) |
| Analysis window hooks: *Project Planner…*, *Script wizard…*, `apply_inputs` (*Use for analysis*: topology / trajectories into Input & cleaning, cleaning on), `run_scripts` (confirmation listing the scripts, background job, log) | app.py:1595-1630, 1944-1987 | - |

MDATK 0.0.11 moved planner and wizard into RNAdo; M5 D1 decided that M7 builds them there.

### 1.10 MDATK tests that specify the workflow layer

| File | Tests | Counterpart |
|---|---|---|
| test_workflows.py | 39: naming, option coercion and requirements, normalise, sequence / dot-bracket parsing, FASTA / secstruct files, `fortran_float`, the Amber stage files against the reference protocol, Amber plan names, build files / errors / warnings, Mg tables and files and checks, GROMACS keeps Mg pairs, Amber MD checks / segments / engines, restraints from a topology and the helper, `prep_pdb`, FARFAR2 files and flags, FARFAR2 models with a fake Rosetta, SimRNA, RNAMake, GROMACS mdp / checks, `amber2gmx` pairs, `gmx_posre`, VMD and groups, connect and order, collisions / cycles / unknown tools, write + master + summary + stamp, every script passes `bash -n` and `py_compile`, plan saved in the project, CLI, Mg pairs and CP types, CP library and TIP4P/2005, CP helper on AmberTools (needs AmberTools), 12-6-4 with CP and two-pass, automatic L-RNA fixes | W2 packages (GoogleTest + generated cases); the AmberTools test as an optional job |
| test_refine.py | 21 (Phenix both versions, defaults, `bash -n`, eff from dot-bracket / model numbering / annotation, checks, L variants need dictionaries, Amber / cpptraj names flagged, builder data, Coot both versions, Python-2 syntax subset, a fake Coot, missing functions, zones, checks, planner inputs, describe / resources, 5'-phosphate flag, inline comments / mmCIF) | W2-D |
| test_workflow_gui.py | 5 (every tool in the stand-alone wizard, stage editor, planner template / connect / write, delete and project round trip, scroll position kept) | W3-B, W3-C (Qt Test; GUI cases) |
| test_workflow_integration.py | 4 (Amber chain runs, Mg pairs reach GROMACS, CP + TIP4P/2005 + nanoMg, 12-6-4 + CP; need AmberTools) | optional job (risk 3) |
| test_projectfile.py | 2 of 13: `test_runner_runs_logs_and_stops`, `..._keyboard_interrupt_and_keeps_output` | W1-C |
| test_cli.py | `test_apply_override` | W2-E |
| test_builder.py, test_dictionaries.py, test_pdbtools.py | the planner-step tests (2, 5 + 1 end to end with a fake acedrg, the embedded copy) | W2-C, W2-D, D3 |

### 1.11 Behaviour confirmed by running MDATK

Probes (scratch scripts against v0.0.4.1, and 0.1.2 where named; Python 3.13, numpy / scipy only;
MDATK's test files `test_workflows.py`, `test_refine.py`, `test_workflow_integration.py`,
`test_cleaning.py`, `test_projectfile.py`, `test_v002.py`, `test_workflow_gui.py`: 94 passed, 9
skipped (no AmberTools, cpptraj or GTK)):

- every template written with stamp `20261127-1200`: steps with errors are skipped, the master keeps
  them as `#   (not written: the step has errors)` and still calls the steps after them (section 5,
  WG-02); `farfar2-amber` with a typed sequence (MDATK's `test_all_generated_scripts_pass_bash_syntax_check`
  settings) has no errors and writes 25 files, master and summary included; `amber-build` with default settings
  writes `cp-lib.py` (149 lines), `prep.py` (107), `pdbtools.py` (1 210), `tleap.leap` (18),
  `parmed.in` (38), `amber-build.sh` (60);
- one-change-at-a-time variants of every generating step type (each check box flipped, each choice
  taken; inputs connected to dummy paths, as 0.1.2's tests/golden.py does): 177 variants, generated
  in 0.1 s in all;
- the `workflow` commands with `--dry-run`, `--set` and `--input` (section 5 lists the failures);
- the two 0.0.4.1 example projects are empty files (`examples/aptamer_target.mdatk.toml`,
  `rna_plan.mdatk.toml`, 0 bytes; restored in 0.0.8, CHANGELOG.md:1887-1888), so there are no example
  plans to compare.

---------------------------------------------------------------------------------------------------

## 2. MDATK after 0.0.4.1 (to 0.1.2): workflow-related changes

Almost all of MDATK's work after 0.0.4.1 went into the workflow layer (the 0.1.1 comparison report,
docs/reports/2026-11-25-mdatk-0.1.1-fix-comparison.md:13): 0.1.2 has 26 step types and 7 templates,
`workflows/` 17 620 lines and `gui/` 12 126 (from 8 684 and about 3 500). Classes: **M7** = adopt now
(a fix of 0.0.4.1 behaviour, a WG item naming the version); **M7-decided** = in M7 by an earlier owner
decision; **later** = the post-parity phase that adds the latest MDATK's features; **convention** =
kept as 0.0.4.1 has it. Parity tests stay on v0.0.4.1, with the M7 fixes applied to MDATK on the test
side (`improvements_workflow.py`).

| Change | Version (CHANGELOG line) | Class |
|---|---|---|
| Hand-picked paths with spaces passed as one argument (build preparation, restraint helper, pmemd `-c` / `-ref`, GROMACS conversion and grompp, SimRNA, RNAMake, FARFAR2 models) | 0.1.0 (CHANGELOG.md:242-244) | M7 (WG-01) |
| `workflow run` reports a step with errors, does not start the steps that read from it, returns 1; Cancel / Stop also end programs a finished script left running | 0.1.0 (:276-278) | M7 (WG-02, WG-16) |
| Temperature ramp of a stage split into segments ran again in every segment; a stage named like another stage's segment or continuation is an error | 0.1.0 (:228-230) | M7 (WG-04, WG-05) |
| Default restraint mask `:1-N` counts residues as tleap reads them | 0.1.0 (:231-233) | M7 (WG-06) |
| FARFAR2 / SimRNA run scripts end with an error when every run failed; SimRNA clustering no longer stops silently (`ls \| head` under pipefail); stale outputs of an earlier run removed first | 0.1.0 (:250-257) | M7 (WG-07 ... WG-09) |
| A secondary-structure file without a sequence is an error that says so | 0.1.0 (:245-246) | M7 (WG-10) |
| Output files written by two steps are an error (`output 'best': ... is also written by ...`) | 0.1.x (plan.py:546 of 0.1.2) | M7 (WG-03) |
| Windows: a save retries while another program holds the file; no `.tmp` left (MDATK-BUG-105) | 0.1.2 (:66-67) | M7 (WG-14, with PRJ-13's helper) |
| Optional number fields keep every digit, also in the script wizard (MDATK-BUG-002) | 0.1.2 (:53-55) | M7 (the M5 form widgets, FW-02) |
| L-RNA and other non-standard atom types get the Mg²⁺ pair terms of their D-RNA twins (*L-RNA atom types*: LD and LE as OS, LH as OH); the microMg / nanoMg frcmod defines both `Mg2+` and `MG`; the Mg block only adds `changeLJPair` lines | 0.0.5 (:1930-1937) | M7 (D9: a physics fix; 0.0.4.1 only warns, amber.py:680-688) |
| *Separate HO and O2 Lennard-Jones types first* on by default for new projects (NBfix19's HO pairs otherwise reach every water hydrogen); two-pass builds leave Mg²⁺ lines out of the dry pass | 0.0.5.1 (:1914-1923) | M7 (D9) |
| Lanthanide trap: an L-RNA library whose units are named otherwise and no `addPdbResMap` makes tleap build LA / LU as ions; the planner reports every unmapped L name with the map lines to add | 0.0.8 (:1872-1876), 0.0.10 (:1751) | M7 for the planner check (D9); the run-time check script later |
| *Check the final topology* (neutral box, Mg²⁺ count and charge, pair terms present; `topcheck.py`) | 0.0.5 (:1926-1929) | later (new option and helper; D9) |
| RNAdo: planner, wizard, PDB clean-up and builder in a program of their own behind a host interface; `rnado [PROJECT] [--open wizard\|pdb\|builder]`; analysis window *Open in RNAdo…* | 0.0.11 (:1630-1668) | M7-decided (M5 D1) |
| Both programs on one project: RNAdo writes only the plan, re-reads before saving, takes newer analysis settings; the analysis window notices and merges | 0.0.11 (:1669-1676, 1697-1746) | M7-decided (M5 plan line 235) |
| Scripts started from a desktop-menu launch get the user's shell set-up (conda) and their own `PYTHONPATH` | 0.0.11, 0.0.13, 0.0.14 (:490-499, :977-986, :1712) | M7 (runner environment, W1-C) |
| Run a selection of steps; independent steps at the same time (*Parallel*), later CPU / GPU budgets and stop buttons | 0.0.6, 0.0.7, 0.0.12 (:1906-1911, :1892-1903, :1367-1410) | later (keys already kept, PRJ-09) |
| Previous runs (*Check previous runs*, `adopted`), completed runs moved aside (`supersede_runs`, `workflow supersede`), output log per script, folder tags, `workflow status / manifest / jobs` | 0.0.10-0.0.13 | later |
| Cluster job scripts, cluster profiles, estimated wall times, `<project>_00_submit_<stamp>.sh`, `mdatk cluster` | 0.1.0 (:167-200) | later |
| Missing programs reported before a step runs; setup file | 0.0.15 (:400-414) | later |
| New step types: `namd-mdff`, `user-script`, `model-check`, `rosetta-audit`, `frame-check`, `register-census`, `ensemble-compare`, `map-model-score`, `body-hybrids`, `decoy-interface`; templates `amber-mdff`, `decoy-validation`; file sets (`collection` ports) and kind `table` | 0.0.8-0.0.15 | later (D10) |
| Amber: generated L-RNA / L-DNA libraries, DNA force fields, RNA force field *None*, per-stage dry trajectories, production stage options, three-character names only | 0.0.8-0.1.0 | later |
| Wizard *Commands* tab; planner layout (boxes keep apart, outlined connections, scrolling while dragging, select all, write selected); cat-proof mode | 0.0.12-0.0.15 | later |
| File naming, step folders, master and summary names, roles, `finished()` stage skipping, option keys, CLI names and JSON | 0.0.4.1 ... 0.1.2 (unchanged; 0.0.8's golden test pins 0.0.7's output) | convention |

---------------------------------------------------------------------------------------------------

## 3. Clean-room boundary, data and licences

### 3.1 MDAnalysis

No workflow module, helper or GUI part imports MDAnalysis (grep, read). The one MDAnalysis-backed path
near the workflow is *Use for analysis*, which only copies paths into the project (app.py:1611-1630).
So milestone 7 needs **no clean-room specs**; implementers may read and run all of MDATK's workflow
code and never open MDAnalysis.

### 3.2 Lua runtime

| Candidate | Licence | Notes |
|---|---|---|
| **Lua 5.4** (PUC-Rio), the newest 5.4 release at least two weeks old when W0-B starts (5.4.7 at the time of writing) | MIT | about 30 C files; can be compiled as C++ (errors then unwind with exceptions, so C++ destructors in the binding run; MDBridge's project is `LANGUAGES CXX`); `lua_sethook` instruction counts, a custom allocator for a memory cap, libraries opened one by one |
| LuaJIT 2.1 | MIT | faster, but Lua 5.1 semantics (no integer subtype, which byte-exact `int` / `float` output needs), rolling releases without version numbers |
| Luau | MIT | sandboxing built in, but a dialect of its own (typed Lua), C++17, large |
| Bindings: sol2 3.3 (MIT, header only) or a hand-written layer | MIT / own | the API is value in, value out (section 4.4), so a small own layer over `core::Value` suffices |

Plan (D1): vendor Lua 5.4 into `third_party/lua` (unmodified sources, pinned by SHA-256 and checked at
configure time, `tools/vendor_lua.py`, as `third_party/ducc0`), compiled as C++ outside the
warnings-as-errors set, licence in `packaging/licenses/lua-LICENSE.txt`, a LICENSING.md row; an own
binding in `src/workflow/lua/`.

### 3.3 Data tables (CLAUDE.md: from `data/sources/` only)

| Table | MDATK's source | M7 |
|---|---|---|
| Force fields (leaprc names), water models (leaprc, box, ParmEd 12-6-4 name, 12-6-4 ion frcmod), four-site waters, default strip masks | amber.py:35-61 ("checked against AmberTools 26 dat/leap") | file names of AmberTools (facts, not content): recorded in `data/sources/amber/leap_files.tsv` with the AmberTools version they were checked against; the generator checks MDATK's tables equal it |
| NBfix19 and 0BPh `changeLJPair` lines and comments | amber.py:63-113 (published parameters: Kührová & Mlýnský et al. 2019; Mlýnský et al. 2023 / 2026) | numbers from the publications in `data/sources/nbfix/`; the comment text is MDATK's own |
| microMg / nanoMg sets (ion, Cl⁻, Mg-Cl and Mg-RNA pairs per water model) | mg_params.py:90-386, read "verbatim from the authors' GROMACS parameter files" at named commits of github.com/bio-phys (Magnesium-FFs, optimizedMgFFs, MgFF_OPC) | `tools/extract_mg_params.py` reads those files at the same commits into `data/sources/mg/*.tsv` (licences of the three repositories to be recorded; the numbers are published in the cited papers either way); `gen_workflow_tables.py` checks MDATK's table equals them (D5) |
| CP values (Steinbrecher, Latzer & Case 2012, Table 1), TIP4P/2005 (Abascal & Vega 2005) | scripts/cp_retype.py:27-29, tip4p2005.py:27-31 | the helpers' text is MDATK's (section 4.6); the numbers checked against the papers in `data/sources/` |
| Hydrogen-bond atom pairs for helix restraints | rna.py:21-30 | chemistry facts, MDATK's own table |
| `FILE_KINDS`, `CATEGORIES`, templates | schema.py:216-234, registry.py:11, plan.py:525-596 | MDATK's own (done for the first two in M5) |

### 3.4 External programs

Generated scripts call tleap, ParmEd (LGPL), ambpdb, pmemd / sander, cpptraj, GROMACS (LGPL), acpype,
Rosetta (`rna_denovo`, `rna_cluster`, `extract_pdbs`; Rosetta licence), x3dna-dssr (own licence),
SimRNA, RNAMake, Phenix, Coot, acedrg, eLBOW, VMD and ffmpeg on the user's machine. MDBridge copies no
code of any of them; their input formats and command lines are behaviour, recorded in the generator
cases. The Coot script (refine.py) and the embedded helpers are MDATK's own code.

---------------------------------------------------------------------------------------------------

## 4. Parity strategy

### 4.1 What is compared

| Artifact | Same means |
|---|---|
| Every generated file (scripts, program inputs, helpers, master, summary) | byte for byte after masking the program name and version in "written by ... / copied by ..." texts and dates (YYYY-MM-DD), as 0.1.2's tests/golden.py normalises (`VERSION_RE`, `DATE_RE`), and the plan root; files a WG item changes: MDATK patched the same way (`improvements_workflow.py`) |
| File names, roles, descriptions, executable bits, the set of files per step | equal |
| Issues (level, message, keys) of every case, `describe`, `resources`, port values (paths, kind, data) | equal |
| Project files after `workflow write` (stamp), `template`, GUI actions | equal text after MDBridge's save rules (PRJ-09, ED-02) |
| CLI output (`workflow tools` text, JSON of `check` / `write` / `generate`, stderr issue lines, exit status) | equal |
| Step-log entries of `workflow write` / `run` | equal entries, time and output masked where they carry time |
| Generated scripts | `bash -n`; helpers `py_compile` (MDATK's test, kept) |

### 4.2 Generators (`tests/parity/`, `regenerate_all.py`)

- `make_workflow_cases.py`: a port of 0.1.2's tests/golden.py case design (owner's code) to 0.0.4.1,
  keeping full texts: every template; every step type with defaults and each check box / choice
  changed one at a time (177 variants, 1.11); each record template of the repeated sections (9 Amber
  stage templates, 5 GROMACS); segments, production and continuation stages; typed secondary
  structures with G·U and U·G (wobble off / on); real input files where checks read them (the PDB
  fixtures with and without L residues for the L-RNA checks, `topology/demo_dry.parm7` for the direct
  DISANG file, models for the dictionary codes and the Phenix eff); planner messages, port values,
  master and summary. Output `tests/data/workflow/cases/*.json` plus the files as text.
- `improvements_workflow.py`: the WG patches applied to MDATK on the test side (as `improvements_*.py`),
  each naming the MDATK version whose design it follows.
  Masking as built (W0-A, accepted): the version after the program's name (`MDBridge <VERSION>`),
  dates (`<DATE>`) and the plan folder (`<ROOT>`) are masked; the name MDBridge itself stays in the
  expected text (D11 sets it), so a file that names the wrong program fails the comparison.
- `make_workflow_goldens.py` (exists, M5): extended with plans that hold every step type's
  non-default settings, saved stamps and scratch records.
- `cli_parity.py`: cases for every `workflow` command and option of 1.8, `--set` records, errors, the
  step log of `write` / `run` (a fake `bash` script set that succeeds, fails and is cancelled).
- `make_gui_cases.py` (M5 harness, GTK under Xvfb): the five test_workflow_gui scenarios as action
  scripts (template, connect, delete, write; wizard per tool; stage editor) with the saved project and
  written files as results.

### 4.3 Tests in C++

GoogleTest suites per area (`tests/workflow/test_plan.cpp`, `test_lua_host.cpp`, `test_amber.cpp`,
`test_gromacs_vmd.cpp`, `test_modelling.cpp`, `test_refine.cpp`, `test_runner.cpp`) read the JSON
cases; each MDATK test of 1.10 gets a counterpart named after it; every generated script passes `bash
-n` on Linux CI. Qt Test for the planner and wizard (W3). Each case runs twice with different Lua
string-hash seeds and must give the same bytes (risk 1).

### 4.4 The Lua API surface (version 1)

A generator file returns one table per step type. C++ calls it; Lua never writes files, starts
programs, reads the clock or the environment.

| Lua sees | From C++ | Notes |
|---|---|---|
| `step{ name, label, category, tag, description, generates, path_settings, inputs, outputs, input_ports(s), output_ports(s), check(s, ctx), output_values(s, ctx), generate(s, ctx), describe(s), resources(s) }` | `workflow::LuaStep`, registered in `core::registerStepType` | built-in steps take their schema from `core` (step_types.inc, D2); a user step declares `sections` in Lua with the same fields as `core::Option` / `Section` |
| `s` | the normalised settings as a read-only table (`core::Value`: Lua integer = Int, float = Float, records as arrays) | |
| `ctx` | `StepContext` (1.2): `project`, `step`, `stamp`, `tag`, `folder`, `title`, `node_id`, `groups`, `preview`, `today` (the date for headers, from C++), `name`, `own`, `from_own`, `has`, `value`, `inp`, `inps`, `data`, `rel`, `existing` (a handle, not a path the script can open), `source_label` | |
| `mdb.issue(level, message, keys)`, `mdb.file(name, content, role, description)`, `mdb.port(paths, kind, data)` | `core::Issue`, `GeneratedFile`, `PortValue` | |
| `mdb.text`: `fortran_float`, `shell_quote`, `extra_lines`, `extra_flags`, `env_lines`, `script_header`, `require_file`, `repr` (Python repr of values), `format(spec, x)` (the Python format specs MDATK uses: `d`, `02d`, `g`, `.Nf`, `.Ne`, `!r`), `utf8` (length and slices in code points) | C++ (`core::pythonRepr`, `formatG` ...) | numbers never go through Lua's `tostring` / `string.format` (risk 1) |
| `mdb.helper(name, constants, call, doc)` | the stored helper texts (4.6) and `embed`'s layout | |
| `mdb.rna`: `parse_sequence`, `parse_dotbracket`, `build_model`, `stems`, `classify_pairs`, `fasta_text`, `secstruct_text`, `HBOND_ATOMS` | `workflow/rna.hpp` (port of rna.py; also used by `workflow restraints`) | |
| `mdb.chem`: the functions of `chem/generator_api.hpp` (parse_ranges, read model, codes, annotation, lw_kind, helix sequence functions, components, SMILES checks) | M6 | file arguments are `ctx.existing` handles |
| `mdb.selections.parse_residue_ranges`, `mdb.mg` (sets, frcmod text, ParmEd lines), `mdb.amber` tables (force fields, waters, NBfix / 0BPh lines), `mdb.citations.has(resource)` | core / workflow tables (3.3) | |
| `mdb.pubchem(id, timeout)` | the program's PubChem transport (M6 D7) | only when `not ctx.preview` and the step's setting asks |

Determinism rules, enforced by a lint step in CI over the bundled `.lua` files: no `pairs` over maps
whose order reaches output (ordered arrays or `mdb.sorted_pairs`), no `tostring` / `..` of floats, no
`os` / `io` / `debug` / `package` (not loaded at all).

### 4.5 User-added `.lua` steps: sandbox and versioning

- Where (owner decision D8): built-in generators are compiled in (`src/workflow/lua/*.lua` as
  resources). User steps are read, without an opt-in, from `<config>/mdbridge/steps/*.lua`, from the
  project's own `steps/*.lua` beside the project file, and from any further folders listed in the
  settings or in the project (`[workflow] step_folders`, relative to the project file, so a project
  and its steps travel together). Order of precedence: project folders, then the settings' folders,
  then the user folder; a step defined in more than one place is taken from the first and a warning
  names every file. The step log records each user step's file path and SHA-256, so a run can be
  reproduced on another machine.
- Sandbox: a fresh `_ENV` per file with only `mdb`, `string`, `table`, `math` (no `math.random`
  seeding from time), `utf8`; `load` of text only; a memory cap (allocator, e.g. 256 MiB) and an
  instruction budget per call (hook), both reported as the step's issue
  `internal error while checking this step: the generator ran out of time / memory`; errors carry the
  file and line.
- Versioning: each file declares `api = 1` and its own `version`; a newer `api` than MDBridge knows is
  refused with a message naming the file; a user step cannot take a built-in name (MDATK's "bad or
  duplicate tool name"); a project that uses a step type that is not available keeps the node and its
  settings (PRJ-09) and reports MDATK's `unknown step type 'x'`.
- Citations: a user step may declare its own references; they join the citation table for that run.

### 4.6 Embedded helper copies (boundary to milestone 8)

MDBridge bundles no Python interpreter today (grep of CMake, packaging and the program). What
generated scripts carry is Python *text*, run by the `python3` of the machine where the simulation
programs run (AmberTools ships one), as in every MDATK version. Plan:

- the six small helpers (`helix.py`, `prep.py`, `cp-lib.py`, `tip4p2005.py`, `convert.py`, `posre.py`)
  are MDATK's own code without known bugs: their texts are extracted from v0.0.4.1 by
  `tools/extract_workflow_helpers.py` into `src/workflow/helpers/` (checked by SHA-256 against
  MDATK's `inspect.getsource` output) and written byte for byte (D4);
- the copies of `pdbtools.py` and `dictionaries.py` are replaced by calls to MDBridge's own `pdb clean`
  and `dict finish` (D3), so the M6 fixes and the CCD data reach the scripts;
- milestone 8 decides whether the remaining helpers become `mdbridge` subcommands.

---------------------------------------------------------------------------------------------------

## 5. MDATK workflow bugs and fragile behaviour to fix (ids WG-xx)

WF-01 ... WF-03 are taken (the density fit in the window, M6 W3-A), so this milestone uses **WG**.
"0.1.2" says whether MDATK's latest `main` still has it and how it was fixed; where it was, MDBridge
follows that design. Each lands as a DIFFERENCES entry (prefix WG).

| Id | MDATK v0.0.4.1 behaviour (verified) | MDBridge | 0.1.2 |
|---|---|---|---|
| WG-01 | Hand-picked inputs with spaces are written unquoted: `-c /tmp/my dir/sys.rst7 -ref /tmp/my dir/sys.rst7` (ran; amber.py:1512 and the other `ctx.inp` uses without `shell_quote`: amber.py:813, 1029; farfar2.py:284, 295; gromacs.py:361, 367; rnamake.py:140; simrna.py:228) | every path in a shell line through `shell_quote` (unchanged text for ordinary paths); a DISANG path with white space (mdin namelist) is a planner error naming it | fixed in 0.1.0: adopt |
| WG-02 | The master script calls the main script of a step whose inputs come from a step with errors (ran: template `amber`, steps 01 / 02 not written, `bash "03_amber-md/..."` still in the master); `workflow run` and the GUI's *Run…* start such steps too (runner.py:94-106) | `run` reports the step with errors, starts none of its readers and exits 1; the master writes `#   (not run: reads from 01 · Amber · system build, which has errors)` for them | `run` fixed in 0.1.0; the master still calls them (plan.py:662-693): adopt + extend |
| WG-03 | Without step folders two steps with one number write the same output files (ran: SimRNA and FARFAR2 models at step 02 both promise `demo_02_model-01_<stamp>.pdb`; only the info "step number 02 is used by ..."), so a reader of one gets the other's models | an error per colliding output, 0.1.2's text (`output 'best': <path> is also written by <step>; give one of them another file name`) | fixed (plan.py:536-548 of 0.1.2): adopt |
| WG-04 | A stage split into segments repeats its temperature ramp: `heat-cont.in` has `temp0=0.0` and `&wt type='TEMP0', istep1=0, istep2=20000, value1=0.0, value2=300.0` again (ran; amber.py:1489-1495, 1524-1527) | the ramp in the first segment only; later segments hold `temp0` | fixed in 0.1.0: adopt |
| WG-05 | A stage named like another stage's segment run (`eq` in 2 segments and a stage `eq-001`) is accepted (ran: no issue); both write `eq-001.out` / `.rst7`, and `finished` then skips the second stage | an error naming both stages | fixed in 0.1.0: adopt |
| WG-06 | The default positional-restraint mask `:1-N` counts residues by (chain, number) (amber.py:228-245): a duplex with blank chain ids and both strands numbered from 1 gives half the solute (ran: two TER-separated strands of 4 residues count 4) | count as tleap reads (a new residue at a name, chain or number change or a TER), 0.1.2's `count_solute_residues` (amber.py:319-341 of 0.1.2) | fixed in 0.1.0: adopt |
| WG-07 | The FARFAR2 and SimRNA run scripts end with status 0 when every run failed: no `set -e`, failures only `echo`ed, the last command an `echo` or `xargs` (ran: generated `farfar2.sh`, `simrna.sh`) | exit 1 with "every run failed" when no output was written; a warning counting failed runs | fixed in 0.1.0: adopt |
| WG-08 | SimRNA clustering: `REF=$(ls runs/r*/*-000001.pdb 2>/dev/null \| head -n 1)` under `set -eo pipefail` (ran: generated script) can end the script without a message when `ls` gets SIGPIPE | no pipeline into `head` (a shell glob loop) | fixed in 0.1.0: adopt |
| WG-09 | Stale outputs: the cluster loop globs `..._thrs*_clust[0-9][0-9].trafl` (ran: generated script), so clusters of an earlier run with another threshold are rebuilt; FARFAR2 models, RNAMake designs and VMD movie frames of an earlier run stay as outputs | each script removes its own earlier outputs first | fixed in 0.1.0: adopt |
| WG-10 | A secondary-structure file with only the dot-bracket line: `internal error while checking this step: string index out of range` (ran; amber.py:155-201, rna.py:215) | "the secondary-structure file has no sequence line: type the sequence" | fixed in 0.1.0: adopt |
| WG-11 | Float settings accept NaN and infinity (schema.py:85-97: comparisons with NaN are false), and generation then fails with a Python text: `workflow generate amber-build --set buffer=nan` -> `error: not enough values to unpack (expected 2, got 1)` (ran; `fortran_float`, base.py:296-314); MDBridge's `core::Option::coerce` reproduces the acceptance today (src/core/step_types.cpp:295-345) | non-finite values are refused in `coerce`: `buffer: expected a finite number, got nan` (on load: the default with that warning, as other bad values) | not fixed (schema.py unchanged) |
| WG-12 | An integer setting given as `1e400` or `inf` ends the CLI with an `OverflowError` traceback (ran: `--set stages.0.maxcyc=1e400`; schema.py:89-90; cli.py:1014 catches no OverflowError) | already MDBridge's: `cannot convert float infinity to integer` as an issue / `error:` line (step_types.cpp:319-327); recorded as a difference | not fixed |
| WG-13 | `workflow run` (and *Run…*) runs the scripts on disk without checking that they are what the current settings write (runner.py:88-128, read): after a settings change without *Write* the old script runs under the new step's title and citations | before running, each step's files are generated and compared with the disk (dates masked); a difference is an error "step 03: its files on disk were written with other settings; write the plan first" (`--force` runs them) | not handled (runner.py, manifest.py) |
| WG-14 | Plan writes use `<path>.tmp` + one `os.replace` (plan.py:471-478, read): on Windows a scanner or indexer holding the file makes the write fail and leaves `.tmp` (MDATK-BUG-105) | PRJ-13's retry helper for every plan file | fixed in 0.1.2 for project files only |
| WG-15 | *New stamp* within the minute of the current stamp gives the same stamp (planner.py:937-939, naming.py:32-33, read), so the new files overwrite the earlier set although the tooltip says "files written before keep their names" | *New stamp* never returns the current stamp (the next free minute, shown) | not fixed |
| WG-16 | Cancel stops the script's process group but not a program the script started in another session that is still running after the script ended (runner.py:56-69, read; 0.1.0's report) | the run ends every process the script started (the M5 process runner's job object / process group, kept until the last child exits) | fixed in 0.1.0: adopt |
| WG-17 | On Windows hand-picked paths are written as `C:\data\x.pdb` into bash scripts that run on Linux / WSL (plan.py:268-278, read; needs Windows to show) | a warning naming the path: "a Windows path in a script for Linux / WSL; copy the file next to the plan or give its Linux path" | not handled |
| WG-18 | The structure-builder script depends on the environment of the writing program (`$APPIMAGE`, frozen or not; buildstep.py:259-268, read), so the same plan writes different files | `mdbridge` by its path when the program is an AppImage, else `mdbridge` on PATH (D6), recorded in the plan summary | not changed |
| WG-19 | The FARFAR2 models script annotates only two-digit models with DSSR: the loop globs `<project>_<NN>_model-[0-9][0-9]_<stamp>.pdb` (farfar2.py:470-472) while the models are numbered `%02d`, so models 100 on are skipped (ran, W0-A: `n_models = 105`, a stand-in `extract_pdbs` and `x3dna-dssr` on PATH: 105 models written, 99 annotated) | the loop globs `model-[0-9]*` (the same pattern as WG-09's removal) | fixed in 0.1.0: adopt |
| WG-20 | SimRNA without clustering promises the output `trajectory` `<project>_<NN>_simrna-all_<stamp>.trafl` (simrna.py:144-149) but only the clustering script writes it (ran, W0-A: `cluster = false`, a stand-in SimRNA: both runs wrote their `.trafl`, the run script ended with status 0, the promised file is missing), so a step reading the trajectory gets nothing | the run script ends with `# all replica trajectories of all runs`, `cat runs/r*/*.trafl > <all>.trafl` and the check that it was written | fixed in 0.1.0: adopt |
| WG-21 | The GROMACS position-restraint helper (`posre.py`, scripts/gmx_posre.py:15-21) does not take the L-DNA names LDA LDC LDG LDT (and their terminal and N variants) for nucleic acid, so an L-DNA strand runs unrestrained (ran, W0-A: `is_nucleic` gives False for LDA, LDA5, LDT3; MDATK 0.0.4.1 itself names no L-DNA, but MDBridge's names do, NA-01) | 0.1.0's `is_nucleic`: an `LD` prefix before the base letter; every name recognised before is recognised as before (the helper's only change, D4) | fixed in 0.1.0: adopt |

Not bugs, kept: `fortran_float`'s 7 significant digits outside [1e-4, 1e7) (base.py:310-314; no
default reaches it), plan order of independent steps by canvas position, the master's "run separately"
lines, `finished()` stage skipping, Windows writing but not running scripts (D12), the default
protocol and every option default except D9's.

---------------------------------------------------------------------------------------------------

## 6. Work breakdown

Code layout: `src/workflow/` (Qt-free library `mdbridge_workflow`, links `core`, `io`, `chem`):
`plan.{hpp,cpp}` (resolve, write, master, summary, templates, scratch), `context.*` (StepContext,
PortValue), `runner.*`, `rna.*`, `tables.inc` (3.3), `helpers/*.py` (4.6), `lua/host.*` (state,
sandbox, `Value` binding), `lua/api.*` (section 4.4), `lua/steps/*.lua` (the 16 step types and shared
Lua modules `amber_common.lua`, `secstruct.lua`), `lua/lint.*`; `src/app/commands/workflow.cpp`;
`src/app/rnado/` (the second program), `src/app/gui/workflow/` (planner, wizard).

### Wave 0 (independent; start together)

| Package | Content | Acceptance |
|---|---|---|
| **W0-A Generators** | section 4.2 generators and `improvements_workflow.py`, `cli_parity.py` workflow cases (skipped by `needs` until the commands exist), the GUI action scripts, regenerate-job steps | `regenerate_all.py --only <steps>` byte-identical twice |
| **W0-B Lua host** | vendored Lua 5.4 (D1), compiled as C++ on GCC 11 and MSVC 2022; state with the sandbox of 4.5, `Value` <-> Lua conversion (integer / float subtypes kept), errors with file and line, memory and instruction limits, the determinism lint | host tests: limits, refused libraries, conversion round trips, the same bytes under two hash seeds |
| **W0-C Tables and helper texts** | `tools/extract_mg_params.py`, `data/sources/mg/`, `data/sources/amber/leap_files.tsv`, `data/sources/nbfix/`, `tools/gen_workflow_tables.py` -> `src/workflow/tables.inc`; `tools/extract_workflow_helpers.py` -> `src/workflow/helpers/`; ATTRIBUTION and LICENSING rows | tables regenerate byte-identically and equal MDATK's; helper texts equal `embed` output |

### Wave 1

| Package | Needs | Content | Acceptance |
|---|---|---|---|
| **W1-A Plan engine** | M5 core | `StepContext`, `PortValue`, resolve (order, inputs, hand-picked files, checks, collisions), write (atomic, modes, WG-14), master, summary, citation block, templates, scratch / single-step plans, `analysis_inputs`, `apply_tool_override`, graph edits; WG-02 (master), WG-03, WG-15 | the plan-level cases of `make_workflow_cases` with stub generators; test_workflows' naming / connect / collisions / write counterparts |
| **W1-B Lua API v1** | W0-B, W0-C | section 4.4: `mdb.*`, `rna.hpp` (port of rna.py), Python formatting helpers, chem binding, the step registry for Lua steps | API tests; `rna` cases (parsing, FASTA, secstruct, pair classes) equal |
| **W1-C Runner** | W1-A, core::process | `plan_jobs`, `run_jobs`, step-log entries, cancel, the Windows message, environment of desktop launches (2), WG-02 (run), WG-13, WG-16 | test_projectfile's two runner counterparts; fake scripts that pass, fail, hang and spawn children |

### Wave 2 (generators in Lua; packages in parallel)

| Package | Needs | Content | Acceptance |
|---|---|---|---|
| **W2-A Amber** | W1-A, W1-B | `amber-build` (prep, CP library, TIP4P/2005, L-RNA through `mdbridge pdb clean` (D3), Mg sets, 12-6-4, ParmEd, two-pass), `amber-restraints` (+ the C++ `workflow restraints`), `amber-md` (stages, segments, production, engines, dry outputs); WG-01, WG-04 ... WG-06, WG-11, D9 (with D9's `split_types` default in `tools/gen_step_types.py`, `src/core/step_types.inc` and the M5 project goldens in one commit: see D9's note) | every Amber case byte for byte; the Amber tests of 1.10 |
| **W2-B GROMACS and VMD** | W1-A, W1-B | `gromacs` (conversion helpers, mdp, posre, mdrun flags), `vmd` (Tcl, groups, movie); WG-09 (VMD frames), WG-21 | their cases; amber2gmx / posre helper tests run with `python3` |
| **W2-C Modelling** | W1-A, W1-B | `farfar2`, `farfar2-models`, `simrna`, `rnamake`, `structure-builder` (D6; non-WC pairings once M6 W2-A lands); WG-07 ... WG-10, WG-18 ... WG-20 | their cases; the fake-Rosetta test; builder planner-step tests |
| **W2-D Refinement and planning** | W1-A, W1-B | `restraint-dicts` (D3), `phenix-rsr`, `coot-rsr`, `note`, `input-file`, `mdatk` | their cases; test_refine counterparts; the dictionary-step tests incl. the fake acedrg run |
| **W2-E Commands** | W1-A, W1-C (grows with W2-A ... D) | `workflow tools / check / write / run / generate / template / restraints` | every `cli_parity.py` workflow case |

### Wave 3

| Package | Needs | Content | Acceptance |
|---|---|---|---|
| **W3-A RNAdo program** | W1-A, M6 W3-B host | `rnado` executable (AppImage action, Start-menu entry, `rnado [PROJECT] [--open wizard\|pdb\|builder]`), the host interface extended (run / stop scripts, use for analysis), the PDB clean-up and builder windows moved, the analysis window's Workflow menu opening RNAdo, co-editing (RNAdo saves only `[workflow]`, re-reads first; the analysis window merges) | GUI cases; both programs on one file keep each other's changes |
| **W3-B Planner window** | W3-A, W2 | canvas, toolbar, naming bar, details pane, write / run with confirmation, *Open folder*, *Use for analysis*, undo for plan edits (M5's edit model) | the planner GUI cases; Qt Test at 100 / 150 / 200 % scale |
| **W3-C Script wizard** | W3-A, W2, M5 forms | forms and record editor from the schema, files page, issues, live preview, write / run | `test_standalone_wizard_every_tool`, `test_wizard_stage_editor` counterparts |
| **W3-D User steps** | W1-B | loading from the user, project and listed folders with their precedence (D8), settings page listing them with errors, an example step and its documentation | a user step appears in *Add step*, writes files, cannot read files or run longer than its budget |
| **W3-E Docs** | all | DIFFERENCES "Workflow" (every WG item: what, why, how tested, version followed), LICENSING (Lua), ATTRIBUTION (Mg sets, NBfix / 0BPh, CP, TIP4P/2005), ARCHITECTURE status, `docs/specs/m7-status.md` | reviewed against the code |

### Critical path and parallelism

`W0-B → W1-B → W2-A → W3-B / W3-C`: the Lua host and API decide every generator, and the Amber family
is the largest (1 607 lines with the most interacting options and the generator cases with real
input files). W1-A is on every path but has no runtime risk; W2-B ... W2-D are independent of each
other and of W2-A once W1-B is done. Three agents in wave 0, three in wave 1, five in wave 2, four or
five in wave 3. The non-WC builder step waits for M6 W2-A only for its non-WC cases.

Packages that depended on the owner decisions (all decided): W0-B (D1), W0-C (D4, D5), W1-B (D2), W1-C (D12), W2-A (D3, D9),
W2-C (D6), W2-D (D3), W3-A (D7), W3-D (D8), all generators (D11), scope (D10).

---------------------------------------------------------------------------------------------------

## 7. Owner decisions (decided 2026-11-28)

The owner accepted every recommendation below except D8, which was changed as written there.

- **D1 Lua runtime**: recommend Lua 5.4 (MIT) vendored unmodified in `third_party/lua`, pinned by
  SHA-256, compiled as C++, with an own value-in / value-out binding (no sol2, LuaJIT or Luau).
- **D2 What is Lua and what is C++**: recommend C++ for the plan engine, writer, master and summary,
  runner, CLI, shared parsers (rna, chem, selections) and tables; Lua for each step type's ports,
  checks, outputs, files, box text and resources; built-in schemas stay generated in `core` from MDATK
  (step_types.inc), user steps declare theirs in Lua.
- **D3 The pdbtools.py / dictionaries.py copies**: recommend replacing them in generated scripts with
  `"$MDBRIDGE" pdb clean ...` / `dict finish ...` (`MDBRIDGE` set at the top of the script to the
  writing program's path, overridable), so the M6 fixes and the CCD data apply; the alternative is
  embedding Python copies patched with the CH fixes, kept in step with the C++ by tests.
- **D4 The six small helpers**: recommend writing MDATK 0.0.4.1's texts byte for byte (extracted and
  checked by a tool), each changed only by a WG item; milestone 8 decides on `mdbridge` subcommands.
- **D5 Mg²⁺ parameter source**: recommend extracting the microMg / nanoMg sets from the authors'
  GitHub files at MDATK's commits into `data/sources/mg/` (licences recorded), checked equal to
  MDATK's table.
- **D6 How scripts start MDBridge** (structure builder; D3): recommend the program's own path when it
  runs as an AppImage, `mdbridge` on PATH otherwise; a non-default `python` setting keeps MDATK's
  `<python> -m mdatk` (the option key stays).
- **D7 GUI host**: recommend confirming M5 D1: planner, wizard, PDB clean-up and builder in the
  `rnado` program with co-editing, opened from the analysis window's Workflow menu.
- **D8 User `.lua` steps** (owner, 2026-11-28: the recommendation of the user folder only was too
  conservative; projects must be shareable between machines with their steps, and the user organizes
  them as they see fit): loaded without an opt-in from the user folder, the project's `steps/` and
  folders listed in the settings or the project, in that precedence from the project down (4.5);
  the sandbox of 4.5 stays (it also keeps the generators deterministic).
- **D9 Physics fixes of 0.0.5-0.0.8 in M7**: recommend adopting the L-RNA Mg²⁺ pair terms of the D
  twins, both Mg type names in the frcmod, *Separate HO and O2 types* on for new projects, the dry
  pass without Mg lines and the planner's lanthanide-trap error; leave *Check the final topology* for
  the post-parity phase.
  Note (W0-A): *Separate HO and O2 types* on by default changes the `split_types` default of the
  amber-build schema, which `tools/gen_step_types.py` takes from MDATK into `src/core/step_types.inc`, and
  with it the M5 project goldens that hold normalised amber-build settings
  (`tests/data/project-file/workflow/golden_*.mdatk.toml`, written by `make_workflow_goldens.py` without
  the improvements) and any other expected data with a normalised amber-build node. The package that
  implements D9 (W2-A) changes them together in one commit: the default in `gen_step_types.py` (on the
  test side `improvements_workflow.py` already sets it), the regenerated `step_types.inc`, and
  `make_workflow_goldens.py` with the D9 default applied, its goldens regenerated.
- **D10 0.1.x features**: recommend none in M7 beyond D7 and D9 (no cluster jobs, budgets, previous
  runs, supersede, folder tags, new step types, `user-script`); they start the post-parity phase.
- **D11 Program name and dates in generated files**: recommend "written by MDBridge <version>" /
  "MDBridge only writes this script" (as CMD-08 and M6 D13), every other byte as MDATK; parity masks
  name, version and date.
- **D12 Running scripts on Windows**: recommend MDATK's behaviour (written, not run, with its
  message) in M7; running through WSL later.
- **D13 Non-finite settings (WG-11)**: recommend refusing NaN and infinity for every float option,
  also when a project is loaded (the default with a warning).
- **D14 Stale scripts (WG-13)**: recommend `run` refusing steps whose files on disk differ from what
  the settings write, with `--force` to run them anyway.

---------------------------------------------------------------------------------------------------

## 8. Risks

1. **Byte parity from Lua**: Python's `repr`, f-string formats and dict order reach the files
   (`embed`'s constants, mdin and mdp values, issue texts); Lua's `tostring`, `string.format` and
   table order differ. Mitigation: all number and value formatting through C++ helpers (4.4), the
   determinism lint, every case run under two hash seeds, 177+ one-change-at-a-time cases.
2. **Size of the port**: about 7 000 lines of generator logic (1 607 in `amber.py` alone) with
   interacting options. Mitigation: the case design of 0.1.2's golden test (each option alone), the
   Amber family split from the others, generators reviewed against MDATK's tests.
3. **External programs cannot run in CI**: the scripts are compared as text and syntax-checked;
   tleap, pmemd, GROMACS, Rosetta never run (MDATK's own integration tests skip without AmberTools).
   Mitigation: an optional CI job with AmberTools from conda-forge running the `amber` template on the
   demo structure (MDATK's test_workflow_integration counterparts); fake programs for the run paths.
4. **Plan-time file reads**: L-RNA detection, the solute count, the direct DISANG file, dictionary
   codes and Phenix pairs depend on input files, so cases need real fixtures. Mitigation: the
   generator copies them into the case folders; the C++ tests read the same files.
5. **Sandbox and security of user steps** (D8): Lua code in a shared project runs when the project
   is opened. The owner judges this an unlikely attack vector and chose project folders for sharing.
   Mitigation: the sandbox (no `io` / `os` / `package` / `debug`, files only through the binding,
   memory and time limits), the step log naming each file with its SHA-256, and a warning for a step
   defined in more than one place.
6. **GUI scope**: RNAdo, planner canvas, wizard and co-editing are large (0.1.2's GUI has grown to
   12 000 lines). Mitigation: 0.0.4.1's features only (D10), M5's forms, jobs and edit model reused,
   co-editing limited to "RNAdo writes only the plan".
7. **Two implementations of chem behaviour in scripts** (D3 alternative, D4): Python text beside C++
   can drift. Mitigation: D3's call to `mdbridge`, and the helper texts frozen with a checksum.
8. **Later MDATK features after parity** add keys and step types; the schema must stay data, and
   0.1.x project files must keep opening (PRJ-09 already keeps their keys).
9. **Windows**: scripts are bash for Linux; Windows-written plans are only checked by the tests on
   Linux (path rules of ARCHITECTURE "Compatibility"; WG-17).
