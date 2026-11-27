# MDATK bug report, addendum 1: further issues still present in 0.1.1

- **Date:** 2026-11-26
- **Follows:** `docs/reports/2026-11-25-mdatk-bug-report.md` (MDATK-BUG-001 … 085). This addendum
  repeats none of its items; where an item here is a variant of an earlier one, it says so.
- **Source:** new findings since 2026-11-25 while building MDBridge (the C++ successor of MDATK): the chem
  layer review (MDBridge `docs/specs/m6-plan.md` section 5, ids `CH-xx`), the app layer review
  (`docs/specs/m5-plan.md` section 4, `AP-xx`) and the `docs/DIFFERENCES.md` entries added since the first
  report (`PRJ-`, `CMD-`, `WA-`, `LE-`, `JB-`, `PKG-` …).
- **MDATK version checked:** the newest `origin/main` after `git fetch` on 2026-11-26:
  commit `1efa0ef44e9cb0d04b9c04a409919f7913f59fdb` ("Claude/exciting euler q7bx2f (#31)", 2026-11-25),
  `pyproject.toml` version **0.1.1**. This is the same commit the first report was checked against; the
  fetch only moved the feature branch `claude/exciting-euler-q7bx2f`, whose tree is identical to `main`
  (`git diff origin/main origin/claude/exciting-euler-q7bx2f` is empty).
- **Test environment:** as the first report: Python 3.13.18, MDAnalysis 2.11.0, numpy 2.5.5, scipy 1.18.1,
  matplotlib 3.11.2; GUI items with PyGObject 3.48.2 / GTK 3 under Xvfb (`xvfb-run -a`). MDATK run from the
  checkout (not installed) with `python -I -B` and the checkout on `sys.path`.

## How to use this addendum

Same conventions as the first report:

- Ids continue at `MDATK-BUG-086`. The MDBridge register id is for cross-reference only.
- `file:line` is MDATK 0.1.1 (the commit above), relative to the repository root.
- **Status:** *Verified* = reproduced against 0.1.1 with the command shown (output quoted as observed);
  *Confirmed by reading* = the cited code still does what is described (not run).
- **Inherited from MDAnalysis:** none of the items below is; all are in MDATK's own code.
- "Suggested fix" says what MDBridge did (or plans, for the chem items of its milestone 6), in words, and
  why. Prefer whatever fits MDATK's design.
- Reproductions assume a shell in the repository root, `mdatk` meaning `python -m mdatk`, and Python
  snippets run with the checkout on `sys.path`. Several use the demo system
  (`mdatk demo demo --n-frames 6 --quiet`) or a small helix and map:
  `mdatk build helix --seq GGCAUGCC -o h.pdb`, and `h.mrc` / `hnan.mrc` from the snippet of 091.
- Severity scale: **critical data loss** > **wrong results** > **crash** > **misleading message** >
  **minor**.

## Summary

21 items: 2 critical data loss, 5 wrong results, 0 crash, 6 misleading message, 8 minor.
17 Verified, 4 Confirmed by reading.

| Id | Area | Severity | Title | Status |
|---|---|---|---|---|
| MDATK-BUG-086 | Project files | critical data loss | Unknown keys (a later MDATK's, another tool's) are deleted by the next save | Verified |
| MDATK-BUG-087 | CLI `density fit --save` | critical data loss | Saving the fit writes every path of the project absolute (variant of 006) | Verified |
| MDATK-BUG-088 | PDB tools / builder | wrong results | Residue numbers below −999 shift the PDB columns; MDATK's own reader then rejects the file | Verified |
| MDATK-BUG-089 | PDB tools | wrong results | `--combine-into XY` writes chain `X` while the notes say `XY` | Verified |
| MDATK-BUG-090 | PDB tools | wrong results | `--keep-conect` drops every CONECT record when input serials repeat, notes say "Kept" | Verified |
| MDATK-BUG-091 | Density fit | wrong results | One NaN voxel is smeared over 729 voxels and changes the fit (variant of 060) | Verified |
| MDATK-BUG-092 | `demo` | wrong results | `--dt` ≤ 0 / NaN writes a trajectory going back in time or at NaN ps; `--n-frames 0` fails after writing | Verified |
| MDATK-BUG-093 | Builder | misleading message | `build helix` wraps residue numbers above 9999 silently | Verified |
| MDATK-BUG-094 | Density | misleading message | Resolution NaN / ∞ and negative spacing end in `negative dimensions are not allowed` | Verified |
| MDATK-BUG-095 | PDB tools / builder | misleading message | A missing input file is reported as "no ATOM/HETATM records found", without its name | Verified |
| MDATK-BUG-096 | PDB tools | misleading message | A reversed range `A:4-1` is reported only as "the L ranges match no nucleotide" | Verified |
| MDATK-BUG-097 | GUI layout page | misleading message | Panel data labels keep the old groups after an analysis' groups change | Verified |
| MDATK-BUG-098 | GUI | misleading message | *Open video* plays the last session render, even of another project or output folder | Verified |
| MDATK-BUG-099 | Builder `mutate` | minor | Mutating a model with hydrogens leaves the new base without hydrogens (mixed model) | Verified |
| MDATK-BUG-100 | Cleaning | minor | Cancel waits for cpptraj's next output and leaves `conda run`'s cpptraj running | Verified |
| MDATK-BUG-101 | Cleaning | minor | cpptraj output decoded per 4 KiB chunk: split UTF-8 characters become two U+FFFD | Verified |
| MDATK-BUG-102 | GUI log | minor | One main-loop callback per log line: a chatty job stalls the analysis window | Verified |
| MDATK-BUG-103 | GUI dialogs | minor | Dialog details still cut at 3 000 characters silently; still preview fixed at 1 200 px | Confirmed by reading |
| MDATK-BUG-104 | GUI | minor | After *Run cleaning now* the preview keeps showing the raw frame | Confirmed by reading |
| MDATK-BUG-105 | Project files (Windows) | minor | A save fails when a scanner or indexer holds the file for a moment; `.tmp` is left behind | Confirmed by reading |
| MDATK-BUG-106 | AppImage | minor | One missing Wayland library puts every bundled fallback library ahead of the host's | Confirmed by reading |

---------------------------------------------------------------------------------------------------

## 1. Critical data loss

### MDATK-BUG-086 — Unknown keys are deleted by the next save

- **MDBridge ref:** PRJ-09 (owner decision D9 of MDBridge milestone 5) · **Area:** project files ·
  **Severity:** critical data loss · **Status:** Verified
- **Affected:** 0.0.4.1 … 0.1.1
- **Where:** `mdatk/config.py:1017-1022` (`_from_plain`: a key the dataclass does not know is logged and
  skipped), `:533` (the meta keys, `mdatk_version` included, are dropped on read);
  `mdatk/projectfile.py:276-279` (`project_text` writes the current `mdatk_version`).
- **Description:** every key MDATK's model does not know is dropped when a project is read, so the next
  save (GUI *Save*, `density fit --save`, and RNAdo's plan saves, which re-read and rewrite the file through
  the same model: `mdatk/projectfile.py:365-378`, `edit_workflow`) deletes it. This loses the settings
  of a newer MDATK (e.g. a file of a future release opened in 0.1.1, as 0.0.4.1 lost 0.1.x's plan budgets,
  adopted runs and folder tags), of other tools that keep their own table in the file, and typos the user
  would have fixed. The warning says "ignoring", not that the key will be removed; in the GUI the warnings
  of the start-up load go to stderr only. The file's `mdatk_version` is replaced by the version that saved.
- **Reproduce:**
  ```python
  import os, tempfile, logging
  logging.basicConfig(level=logging.WARNING, format="%(levelname)s: %(message)s")
  from mdatk.config import Project
  d = tempfile.mkdtemp(); src = os.path.join(d, "p.mdatk.toml"); Project().save(src)
  text = open(src).read().replace("[inputs]\n", "[inputs]\nfuture_option = 3\n", 1)
  open(src, "w").write(text + "\n[plugin_x]\nkeep_me = \"yes\"\n")
  p = Project.load(src); p.save(src); out = open(src).read()
  print("future_option kept:", "future_option" in out, "| [plugin_x] kept:", "plugin_x" in out)
  ```
- **Actual:**
  ```
  WARNING: project.inputs: ignoring unknown key 'future_option'
  WARNING: project: ignoring unknown key 'plugin_x'
  future_option kept: False | [plugin_x] kept: False
  ```
- **Expected:** a save keeps what the program did not understand.
- **Suggested fix:** MDBridge keeps every unknown key of every table at any depth (lists of tables
  included) with its value and type as read, writes it back after the table's own fields, and warns once:
  `project: 2 unknown keys kept and saved back unchanged (from a newer MDATK, or misspelt?): project.bogus,
  project.inputs.nope`. `--set` still refuses keys the model does not know. `mdatk_version` stays what a
  later version wrote. In MDATK: stash the unknown items per dataclass instance (e.g. a `_extra` dict
  filled in `_from_plain`) and merge them back in `_to_plain`; leave `mdatk_version` alone when it is
  newer than `__version__`.
- **Test idea:** load + save a file with unknown keys at the top level, inside `[inputs]`, inside a
  `[[groups]]` entry and inside an analysis section; assert each key, value and type survives.

### MDATK-BUG-087 — `density fit --save` writes every path of the project absolute

- **MDBridge ref:** CH-12 · **Area:** command line (`density fit`) · **Severity:** critical data loss ·
  **Status:** Verified
- **Affected:** 0.0.4.1 … 0.1.1. Variant of **MDATK-BUG-006** (the GUI save), with its own code path.
- **Where:** `mdatk/cli.py:165-169` (`_pipeline` → `load_project(args)`, which resolves paths by default,
  `:97`), `:450-455` (`_density_overrides` stores `os.path.abspath(args.map)`), `:477-483`
  (`project.save(target)` saves the resolved project).
- **Reproduce:**
  ```
  mdatk demo demo --n-frames 6 --quiet
  mdatk density simulate -p demo/demo.mdatk.toml --resolution 8 --output demo/sim.mrc --quiet
  mdatk density fit -p demo/demo.mdatk.toml --map demo/sim.mrc --resolution 8 --save --quiet
  grep -n '^topology\|^map_path\|^out_dir' demo/demo.mdatk.toml
  ```
- **Actual:** before: `topology = "demo.parm7"`, `trajectories = ["demo.nc"]`, `out_dir = "MDATK_output"`;
  after: `topology = "/…/demo/demo.parm7"`, the trajectories, `map_path = "/…/demo/sim.mrc"` and
  `out_dir = "/…/demo/MDATK_output"` all absolute.
- **Expected:** paths as the file had them; a new `--map` given relative is stored relative to the project
  folder. A project folder that is moved or shared breaks after one `--save`.
- **Suggested fix:** MDBridge saves through its edit model, which keeps every path as written and stores a
  new path given relative to the current directory relative to the project folder (the planned fix for
  CH-12; the GUI part, ED-02, has landed). In MDATK: load with `resolve=False` for the copy that is saved
  (as `init` does, `cli.py:228`), apply the fit result to it, and save that; or the general fix proposed
  for 006.
- **Test idea:** the commands above; assert the three keys' text is unchanged and `map_path` is
  `sim.mrc`.

## 2. Wrong results

### MDATK-BUG-088 — Residue numbers below −999 shift the PDB columns

- **MDBridge ref:** CH-01 · **Area:** PDB tools, builder · **Severity:** wrong results · **Status:** Verified
- **Affected:** 0.0.4.1 … 0.1.1
- **Where:** `mdatk/pdbtools.py:640-643` (`_id_columns`: `f"{resseq:>4}"` widens to 5 characters for
  −1000 and below; numbers above 9999 wrap, numbers below −999 do not). Reached by
  `build helix --start=-1002,1`, `pdb clean --renumber --start -2000`, `--copy-numbers` from such a file,
  `renumber_residues(start=-2000)`.
- **Reproduce:**
  ```
  mdatk build helix --seq GCAU --start=-1002,1 -o n.pdb
  grep " C1'" n.pdb | head -1
  mdatk pdb clean n.pdb
  ```
- **Actual:** `ATOM      9  C1'   G A-1002       6.837   5.134   2.445` (every later column moved one to the
  right), and MDATK's own reader then fails: `error: line 5: unreadable coordinates`.
- **Expected:** no file that MDATK (or any PDB reader) cannot read is written.
- **Suggested fix:** MDBridge's plan: refuse before anything is written, with
  `residue number -1002 (A) does not fit the PDB columns (-999 to 9999)`; keep MDATK's wrap above 9999 (a
  convention) with its warning (see 093). In MDATK: raise `PDBError` in `_id_columns` (or once in
  `write_pdb`) for `resseq < -999`.
- **Test idea:** `build_helix("GCAU", start=(-1002, 1))` raises `PDBError`; `start=(-999, 1)` writes a file
  that `read_pdb` reads back with the same numbers.

### MDATK-BUG-089 — `--combine-into XY` writes chain `X` while the notes say `XY`

- **MDBridge ref:** CH-03 · **Area:** PDB tools (`pdb clean`) · **Severity:** wrong results · **Status:** Verified
- **Affected:** 0.0.4.1 … 0.1.1
- **Where:** `mdatk/pdbtools.py:1659-1680` (`combine_chains` takes `into` as given), `:643` (`[:1]` cuts the
  chain id when writing); option `:2078`. `restore_chains` already refuses such ids (`:1594-1595`,
  `chain ids must be single characters`).
- **Reproduce:**
  ```
  mdatk build helix --seq GCAU -o w.pdb
  mdatk pdb clean w.pdb --combine A,B --combine-into XY -o c.pdb
  grep '^ATOM' c.pdb | cut -c22 | sort | uniq -c
  ```
- **Actual:** notes `Combined chains A, B into chain XY (8 residues)` and `WARNING: chain XY now has duplicate
  residue numbers …`; every ATOM line has chain `X` (164 of 164).
- **Expected:** an error for a chain id that does not fit column 22.
- **Suggested fix:** MDBridge's plan: `PdbError("chain ids must be single characters")`, the message
  `restore_chains` gives, checked before any change.
- **Test idea:** `clean(text, combine=["A", "B"], combine_into="XY")` raises; `"X"` works.

### MDATK-BUG-090 — `--keep-conect` drops every CONECT record when input serials repeat

- **MDBridge ref:** CH-04 · **Area:** PDB tools (`pdb clean`) · **Severity:** wrong results · **Status:** Verified
- **Affected:** 0.0.4.1 … 0.1.1
- **Where:** `mdatk/pdbtools.py:702-704` (`_remap_conect` returns `[]` when any input serial maps to more
  than one atom), `:2013-2015` (the note says `Kept CONECT records (remapped to the new atom serials)`
  whenever there were CONECT records and `strip_conect` is off).
- **Description:** files concatenated from two PDB files, or with wrapped serials, have repeated serials;
  then no CONECT record is written at all, while the notes claim they were kept and remapped.
- **Reproduce:**
  ```
  mdatk build helix --seq GCAU -o w.pdb
  (grep '^ATOM' w.pdb; echo TER; grep '^ATOM' w.pdb | sed 's/^\(.\{21\}\)./\1C/'; echo 'CONECT    1    2'; echo END) > dup.pdb
  mdatk pdb clean dup.pdb --keep-conect -o d.pdb; grep -c CONECT d.pdb
  ```
- **Actual:** `pdbtools: Kept CONECT records (remapped to the new atom serials)`; `0` CONECT lines.
- **Expected:** the note says what happened.
- **Suggested fix:** MDBridge's plan: the note becomes `WARNING: dropped N CONECT records: the input's atom
  serials repeat, so they cannot be remapped`. In MDATK: let `write_pdb` report how many records were
  written (or check the serial map in `clean`) and choose the note from that.
- **Test idea:** the reproduction; assert the warning and that a file without repeated serials still keeps
  and remaps its records.

### MDATK-BUG-091 — One NaN voxel is smeared over 729 voxels and changes the fit

- **MDBridge ref:** CH-11 · **Area:** density fit (`build fit`, `density fit`) · **Severity:** wrong results ·
  **Status:** Verified
- **Affected:** 0.0.4.1 … 0.1.1. Variant of **MDATK-BUG-060** (the density shadow) in the fit.
- **Where:** `mdatk/density.py:340-347` (`FitProblem.at_resolution`: `ndimage.gaussian_filter` of the raw
  data), `:495-503` (the coarse search's anti-alias smoothing, then samples that are NaN count as outside
  the map).
- **Description:** maps from some pipelines carry NaN voxels (masked regions, failed sharpening). The
  Gaussian filter spreads one NaN over its whole kernel (9³ voxels at σ = 1 voxel); every sample there is
  then treated as outside the map, which lowers the inside fraction, the CC and the pose.
- **Reproduce:**
  ```python
  import numpy as np
  from mdatk import builder, density, pdbtools as pt
  from mdatk.io.mrc import write_map
  from scipy import ndimage
  open("h.pdb", "w").write(builder.build_helix("GGCAUGCC"))
  xyz, el = builder.model_atoms(pt.read_pdb(open("h.pdb").read()))
  m = density.simulate_map(xyz + np.array([3.0, -2.0, 1.0]), density.atom_weights(el), 6.0)
  write_map("h.mrc", m)
  d = m.data.copy(); d[tuple(s // 2 for s in d.shape)] = np.nan
  write_map("hnan.mrc", density.DensityMap(d, m.origin, m.voxel))
  print("NaN after sigma=1:", int(np.isnan(ndimage.gaussian_filter(d.astype(np.float32), 1.0)).sum()))
  for mp in ("h.mrc", "hnan.mrc"):
      for search in ("local", "global"):
          _t, cc, notes = builder.fit_to_map("h.pdb", mp, 6.0, search=search)
          print(mp, search, f"CC {cc:.4f}", [n for n in notes if "moved" in n][0])
  ```
- **Actual:** `NaN after sigma=1: 729`;
  ```
  h.mrc    local  CC 0.7864  Model moved by a 29.1° rotation and a 2.9 Å shift of its centre
  h.mrc    global CC 0.9763  Model moved by a 0.0° rotation and a 3.7 Å shift of its centre
  hnan.mrc local  CC 0.7860  Model moved by a 29.1° rotation and a 2.9 Å shift of its centre
  hnan.mrc global CC 0.9763  Model moved by a 0.1° rotation and a 3.7 Å shift of its centre
  ```
- **Expected:** one bad voxel (of 10 920) changes nothing measurable; the user is told about it.
- **Suggested fix:** MDBridge's rule for maps (RN-30 / RN-33, planned for the fit as CH-11): non-finite
  voxels carry no density (0 before any smoothing or resampling, and a mask so they are not counted as
  "inside"), and a note counts them. In MDATK: `np.nan_to_num(data, nan=0, posinf=0, neginf=0)` once
  when the map is loaded for fitting (`load_map_cached`), plus a note with the count.
- **Test idea:** the snippet; assert identical CC and angles for `h.mrc` and `hnan.mrc`, and the note.

### MDATK-BUG-092 — `demo` writes impossible trajectories, or fails after writing

- **MDBridge ref:** CMD-01 · **Area:** `demo` · **Severity:** wrong results · **Status:** Verified
- **Affected:** 0.0.4.1 … 0.1.1
- **Where:** `mdatk/cli.py:1236-1237` (`--n-frames`, `--dt` unchecked), `:202-205`; `mdatk/demo.py:554-568`
  (`build_demo` writes `demo.parm7` before the frames are built), `:202`, `:328` (`np.zeros((n_frames, …))`),
  `:517-524` (`write_trajectory(…, dt=dt_ps)`).
- **Reproduce:**
  ```
  mdatk demo d0 --n-frames 0 --quiet;  ls d0
  mdatk demo dn --n-frames -2 --quiet; ls dn
  mdatk demo dd --n-frames 3 --dt -5 --quiet
  python -c "from scipy.io import netcdf_file as n; print(n('dd/demo.nc','r',mmap=False).variables['time'][:])"
  mdatk demo dnan --n-frames 2 --dt nan --quiet; grep -n nan dnan/demo.mdatk.toml
  ```
- **Actual:** `error: [Errno 5] trying to go over trajectory limit` with `d0/demo.parm7` left behind;
  `error: negative dimensions are not allowed` (empty `dn/`); `--dt -5` succeeds and writes times
  `[-0. -5. -10.]`; `--dt nan` succeeds and the project holds `times = "nanps"`.
- **Expected:** impossible arguments refused before anything is written.
- **Suggested fix:** MDBridge refuses them first: `the demo needs at least one frame (--n-frames 0)` and
  `the time between frames must be a positive number of ps (--dt 0.0)`; the library's sequence must have at
  least four of A, C, G, U. In MDATK: check in `cmd_demo` (or `build_demo`) with `math.isfinite(dt) and
  dt > 0` and `n_frames >= 1`.
- **Test idea:** each command above exits 1 with the message and leaves no files.

## 3. Misleading messages

### MDATK-BUG-093 — `build helix` wraps residue numbers above 9999 silently

- **MDBridge ref:** CH-02 · **Area:** builder · **Severity:** misleading message · **Status:** Verified
- **Affected:** 0.0.4.1 … 0.1.1 (also `build mutate` of such files)
- **Where:** `mdatk/pdbtools.py:640-643` (the wrap); the warning exists only in `renumber_residues`
  (`:1841-1842`), which `build_helix` (`mdatk/builder.py:722`) does not call.
- **Reproduce:** `mdatk build helix --seq GCAU --start=1,9998 -o w.pdb; grep " C1'" w.pdb`
- **Actual:** chain B is numbered `B9998`, `B9999`, `B   0`, `B   1`; the notes say nothing.
- **Expected:** the wrap (a convention MDATK keeps) comes with `WARNING: residue numbers above 9999 wrap in
  the PDB columns`, as `pdb clean --renumber` prints.
- **Suggested fix:** MDBridge adds the same warning to the build notes. In MDATK: emit it from `write_pdb`'s
  callers whenever a number above 9999 is written (or from `build_helix` after numbering).
- **Test idea:** `build_helix("GCAU", start=(1, 9998), notes=n)`; assert the warning is in `n`.

### MDATK-BUG-094 — Resolution NaN / ∞ and negative spacing end in `negative dimensions are not allowed`

- **MDBridge ref:** CH-10 · **Area:** density (`build fit`, `density fit`, `density simulate`) ·
  **Severity:** misleading message · **Status:** Verified
- **Affected:** 0.0.4.1 … 0.1.1. Related to **MDATK-BUG-012** (NaN accepted by validation): the project
  check `mdatk/config.py:673` (`if d.resolution <= 0`) lets NaN through too.
- **Where:** `mdatk/builder.py:2440-2441` (`not resolution or resolution <= 0`: NaN passes),
  `mdatk/density.py:182-189` (`simulate_map`: same test; a negative `spacing` is used as given and the grid
  size becomes negative), `mdatk/cli.py:499` (`args.spacing or None`), `:1319-1320`.
- **Reproduce:**
  ```
  mdatk build fit h.pdb h.mrc --resolution nan -o f.pdb
  mdatk build fit h.pdb h.mrc --resolution inf -o f.pdb
  mdatk density fit -p demo/demo.mdatk.toml --map demo/sim.mrc --resolution nan --quiet
  mdatk density simulate -p demo/demo.mdatk.toml --resolution 6 --spacing -1 --output s.mrc --quiet
  ```
- **Actual:** each prints numpy's `RuntimeWarning: invalid value encountered in cast` (not for the negative
  spacing) and `error: negative dimensions are not allowed`. (`--resolution -3` gives the good message
  `Structure view › Density map › Resolution: Density map resolution must be positive.`; `--step nan` gives
  `the angular step must be at least 3 degrees (got nan)`.)
- **Expected:** "the resolution must be a positive number of Å" / "the grid spacing must be positive"
  before any work.
- **Suggested fix:** MDBridge's plan: those two messages; `--spacing 0` keeps MDATK's meaning (the default
  resolution / 3). In MDATK: `if not (math.isfinite(resolution) and resolution > 0)` in `simulate_map`,
  `fit_to_map` and `density_problems`; `spacing` likewise when given.
- **Test idea:** the four commands exit 1 with the located message and no numpy warning.

### MDATK-BUG-095 — A missing input file is reported as "no ATOM/HETATM records found"

- **MDBridge ref:** CH-13 · **Area:** PDB tools, builder · **Severity:** misleading message · **Status:** Verified
- **Affected:** 0.0.4.1 … 0.1.1
- **Where:** `mdatk/pdbtools.py:407-414` (`_read_text`: a string that is not an existing file is taken as
  PDB text), `:596-597` (the message).
- **Reproduce:** `mdatk pdb clean nosuch.pdb`; `mdatk build mutate nosuch.pdb A:2G -o x.pdb`;
  `mdatk build fit nosuch.pdb h.mrc --resolution 6 -o f.pdb`
- **Actual:** each: `error: no ATOM/HETATM records found (not a PDB file, or the path does not exist)`.
- **Expected:** the command line names the file: `nosuch.pdb: no such file`. (A typo in a long path is hard
  to spot from the current text; a directory or unreadable file gives the same message.)
- **Suggested fix:** MDBridge's plan: the commands open the file themselves and report the missing file
  (its usual file error); the API keeps accepting text. In MDATK: in the command handlers, check
  `os.path.isfile(args.input)` and raise with the path before calling the PDB tools.
- **Test idea:** the commands above print the path and exit 1.

### MDATK-BUG-096 — A reversed range `A:4-1` is reported only as "the L ranges match no nucleotide"

- **MDBridge ref:** CH-15 · **Area:** PDB tools (`--l-ranges`, expected-hand ranges) ·
  **Severity:** misleading message · **Status:** Verified
- **Affected:** 0.0.4.1 … 0.1.1
- **Where:** `mdatk/pdbtools.py:1077-1101` (`parse_ranges` accepts `start > end`), `:1104-1107`
  (`_in_ranges` then matches nothing), `:1181-1182` (and `:1303-1304` for expected-hand ranges).
- **Reproduce:** `mdatk build helix --seq GCAU -o w.pdb; mdatk pdb clean w.pdb --l-ranges A:4-1 -o r.pdb`
- **Actual:** `WARNING: the L ranges match no nucleotide`, exit 0.
- **Expected:** `residue range A:4-1 ends before it starts`, as an error.
- **Suggested fix:** MDBridge's plan: that error from range parsing. In MDATK: raise `PDBError` in
  `parse_ranges` when `b < a`.
- **Test idea:** `parse_ranges("A:4-1")` raises; `"A:1-4"` and `"A:-4--1"` parse.

### MDATK-BUG-097 — Panel data labels keep the old groups after an analysis' groups change

- **MDBridge ref:** LE-06 · **Area:** GUI (Layout & panels) · **Severity:** misleading message · **Status:** Verified
- **Affected:** 0.0.4.1 … 0.1.1
- **Where:** `mdatk/gui/app.py:1641-1651` (`_fill_data` writes `key_label(key, project)` into the slot's
  Data drop-down), only called through `_rebuild_slots` (`:1710`), which runs on load, preset change and
  group-table edits (`_groups_changed`, `:813-825`), not when an analysis' group combo changes
  (`_group_combo`, `:827-830`; contacts at `:1090-1091`).
- **Reproduce** (Xvfb):
  ```python
  import os, tempfile, gi
  gi.require_version("Gtk", "3.0")
  from gi.repository import Gtk
  from mdatk.config import Project, PanelSpec, LayoutParams, ResidueGroup
  from mdatk.gui.app import MainWindow
  app = Gtk.Application(application_id=None); app.register(None)
  src = os.path.join(tempfile.mkdtemp(), "p.mdatk.toml")
  p = Project(); p.groups.append(ResidueGroup(name="Loop", residues="5-8"))
  p.layout = LayoutParams(preset="single", panels=[PanelSpec("timeseries", ["contacts/atom_pairs"], "")])
  p.save(src); win = MainWindow(app, Project.load(src), src)
  pump = lambda: [Gtk.main_iteration_do(False) for _ in iter(Gtk.events_pending, False)]
  combo = win.slot_rows[0]["data"]; label = lambda: combo.get_model()[combo.get_active_iter()][0]
  print("before:", label())
  win.group_combos[1].get_child().set_text("Loop"); pump()   # analyses.contacts.group_b
  print("group_b:", win.collect(validate=False).analyses.contacts.group_b, "| after:", label())
  ```
- **Actual:** `before: Aptamer–Target heavy-atom contacts   [contacts/atom_pairs]`;
  `group_b: Loop | after: Aptamer–Target heavy-atom contacts   [contacts/atom_pairs]`.
- **Expected:** `Aptamer–Loop heavy-atom contacts   [contacts/atom_pairs]`.
- **Suggested fix:** MDBridge's editors refresh their labels after every change (selections and titles
  stay). In MDATK: connect the analysis group combos to a light refresh that relabels the slot drop-downs
  (keeping the active ids) without rebuilding the editors.
- **Test idea:** the snippet; assert the label names `Loop`.

### MDATK-BUG-098 — *Open video* plays the last session render, even of another project or folder

- **MDBridge ref:** WA-03 · **Area:** GUI · **Severity:** misleading message · **Status:** Verified
- **Affected:** 0.0.4.1 … 0.1.1
- **Where:** `mdatk/gui/app.py:2364-2370` (`open_video` prefers `self.last_video` whenever the file still
  exists), `:2397-2399` (set by *Render* only), `:2513-2534` (*Re-render from cache* does not set it),
  `:109` (never reset: `load_project`, *New*, *Close* keep it); `mdatk/gui/desktop.py:97-102`
  (`latest_video`, ties broken by directory order).
- **Description:** after a *Render*, *Open video* keeps opening that file after the output folder was
  changed, after a *Re-render from cache* wrote a newer video, and after another project was opened.
- **Reproduce** (Xvfb; continuing the window of 097):
  ```python
  win.last_video = os.path.join(os.path.dirname(src), "old_project_video.mp4"); open(win.last_video, "w").close()
  src2 = os.path.join(os.path.dirname(src), "other", "q.mdatk.toml"); os.makedirs(os.path.dirname(src2))
  Project().save(src2); win.load_project(Project.load(src2), src2); pump()
  print(win.last_video)
  ```
- **Actual:** the first project's video path; `open_video()` would open it for the second project.
- **Expected:** the newest video of the output folder the current settings name.
- **Suggested fix:** MDBridge remembers the video of a *Render* or a *Re-render* and opens it only while it
  lies in the folder the settings name; otherwise the newest matching file (equal times: the last name;
  only files, not folders named like a video). In MDATK: set `last_video` in `rerender_from_cache`, reset it
  in `load_project` / `new_project`, and check `os.path.dirname(last_video) == project.out_dir`.
- **Test idea:** render into A, change `out_dir` to B containing a video, assert B's video is opened.

## 4. Minor

### MDATK-BUG-099 — Mutating a model with hydrogens leaves the new base without hydrogens

- **MDBridge ref:** CH-09 · **Area:** builder (`build mutate`) · **Severity:** minor · **Status:** Verified
- **Affected:** 0.0.4.1 … 0.1.1
- **Where:** `mdatk/builder.py:2188-2255` (`_swap_base`: the old base's hydrogens are dropped, the new base
  comes from `template(base, hydrogens=False)`).
- **Description:** the note says so (`the new base has heavy atoms only`), but the result is a mixed model:
  every other residue (the sugar of the mutated one included) has hydrogens, the new base has none. tleap
  then adds hydrogens to that base only, and the model checks and H-bond analyses treat it inconsistently.
- **Reproduce:**
  ```python
  from mdatk import builder, pdbtools as pt
  t, notes = builder.mutate(builder.build_helix("GCAU", hydrogens=True), "A:2G")
  r = [r for r in pt.read_pdb(t).residues() if (r.chain, r.resseq) == ("A", 2)][0]
  print(notes[1]); print(r.name, [a.name for a in r.atoms if a.name.startswith("H") and "'" not in a.name])
  ```
- **Actual:** `A:2: dropped 4 hydrogen(s) of the old base (H41 H42 H5 H6); the new base has heavy atoms only`;
  `G []` (7 sugar hydrogens kept; residue A:1, a G, has `H8 H1 H21 H22`).
- **Expected:** a model with hydrogens stays one.
- **Suggested fix:** MDBridge's plan: when the old base had hydrogens, the new base gets them from the
  builder's ideal rules (`_hydrogens`, which `build helix --hydrogens` already uses), with a note.
- **Test idea:** mutate each base into each other one in a `hydrogens=True` helix; assert the new base has
  the hydrogens `template(base, hydrogens=True)` lists.

### MDATK-BUG-100 — Cancel waits for cpptraj's next output and leaves `conda run`'s child running

- **MDBridge ref:** CMD-03 · **Area:** cleaning · **Severity:** minor · **Status:** Verified
- **Affected:** 0.0.4.1 … 0.1.1
- **Where:** `mdatk/cleaning/cpptraj.py:204-216` (the cancel flag is checked only between reads;
  `read1` at `:217` blocks until cpptraj prints; `proc.terminate()` ends the started program, not its
  children; `Popen` without a new session / process group).
- **Description:** cpptraj is silent for long stretches of a large trajectory, so *Cancel* waits until it
  next prints. With `conda run`, terminate ends `conda` and leaves cpptraj running, writing into the
  cleaned folder.
- **Reproduce** (fake programs; `fake_wrapper.sh` stands in for `conda run … cpptraj`):
  ```sh
  printf '#!/bin/sh\necho started\nsleep 2; echo "progress 10%%"\nsleep 23\n' > fake_wrapper.sh; chmod +x fake_wrapper.sh
  ```
  ```python
  import os, threading, time, subprocess
  from mdatk.cleaning.cpptraj import run_cpptraj, Cancelled
  ev = threading.Event(); threading.Timer(1.0, ev.set).start(); t = time.time()
  try: run_cpptraj([os.path.abspath("fake_wrapper.sh")], "x.in", "w.log", ".", lambda l: None, None, ev)
  except Cancelled: print("returned after %.1f s" % (time.time() - t))
  time.sleep(0.3); print(subprocess.run(["pgrep", "-af", "^sleep 23$"], capture_output=True, text=True).stdout)
  ```
- **Actual:** cancel requested at 1.0 s, `returned after 2.0 s` (a wrapper silent for 6 s: 6.0 s); then
  `2088 sleep 23` is still running.
- **Expected:** cancel ends cpptraj (and what `conda run` started) within a fraction of a second.
- **Suggested fix:** MDBridge's process runner polls every 50 ms, starts the program in its own process
  group (Windows: a job object) and ends the whole group: SIGTERM, then SIGKILL after 10 s. In MDATK:
  `Popen(..., start_new_session=True)` (Windows: `CREATE_NEW_PROCESS_GROUP` or a job object), read stdout
  in a helper thread so the loop can poll the event, and `os.killpg(proc.pid, SIGTERM)`.
- **Test idea:** the reproduction; assert return within 0.5 s and no `sleep 23` left.

### MDATK-BUG-101 — cpptraj output decoded per chunk: split UTF-8 characters become U+FFFD

- **MDBridge ref:** CMD-04 · **Area:** cleaning · **Severity:** minor · **Status:** Verified
- **Affected:** 0.0.4.1 … 0.1.1
- **Where:** `mdatk/cleaning/cpptraj.py:217-221` (each 4 KiB chunk is decoded on its own with
  `errors="replace"`, and `cpptraj.log` gets the decoded text).
- **Reproduce:** a fake cpptraj `fake_utf8.py` that writes `b"a" * 4095 + "é\n".encode()`, run through
  `run_cpptraj([sys.executable, "fake_utf8.py"], "x.in", "u.log", ".", lines.append)`.
- **Actual:** the line holds two U+FFFD and no `é` (log file likewise).
- **Expected:** `é` (file names, user paths and titles in cpptraj's output are UTF-8).
- **Suggested fix:** MDBridge decodes whole lines (invalid bytes still become U+FFFD) and writes the log
  with the bytes cpptraj wrote. In MDATK: `codecs.getincrementaldecoder("utf-8")(errors="replace")` across
  chunks, or `io.TextIOWrapper(proc.stdout, encoding="utf-8", errors="replace")` read line by line.
- **Test idea:** the fake program above; assert `é` in the logged line and in `cpptraj.log`.

### MDATK-BUG-102 — One main-loop callback per log line: a chatty job stalls the analysis window

- **MDBridge ref:** JB-02 · **Area:** GUI (analysis window log) · **Severity:** minor · **Status:** Verified
- **Affected:** 0.0.4.1 … 0.1.1 (RNAdo's `LogPanel` of 0.1.x batches; the analysis window does not)
- **Where:** `mdatk/gui/worker.py:14-21` (`call_soon`: one `GLib.idle_add` per call), `:31-33` (each log
  line of a job), `mdatk/gui/app.py:89-97` (`GuiLogHandler`), `:2556-2563` (`append_log` inserts and trims
  per line).
- **Reproduce** (Xvfb; the window of 097):
  ```python
  import threading, time
  from gi.repository import GLib
  from mdatk.gui.worker import call_soon
  th = threading.Thread(target=lambda: [call_soon(win.append_log, f"line {i}") for i in range(200000)])
  th.start(); th.join(); fired = {}
  GLib.timeout_add(10, lambda: fired.setdefault("t", time.time()) and False); t1 = time.time()
  while "t" not in fired: Gtk.main_iteration_do(True)
  print("10 ms timer ran after %.1f s" % (fired["t"] - t1))
  ```
- **Actual:** 200 000 lines queued in 1.0 s; a 10 ms timer scheduled after them ran after 1.5 s (the
  window does not redraw or react meanwhile; it grows linearly with the lines, and a producer that keeps
  printing keeps the loop busy).
- **Expected:** the window stays responsive whatever a job prints.
- **Suggested fix:** use 0.1.x's own design (`mdatk/gui/logpanel.py`: lines queued from any thread, inserted
  in batches every 40 ms, at most 2 000 per batch with a note about skipped lines) for the analysis window
  too; MDBridge does this and coalesces progress updates (at most one pending, the newest value).
- **Test idea:** a job that logs 10⁶ lines; assert no main-loop stall over 0.2 s and that Cancel works
  during the flood.

### MDATK-BUG-103 — Dialog details still cut at 3 000 characters silently; still preview fixed at 1 200 px

- **MDBridge ref:** AP-18, WA-07 · **Area:** GUI dialogs · **Severity:** minor · **Status:** Confirmed by reading
- **Affected:** 0.0.4.1 … 0.1.1. The first report listed AP-18 as fixed in 0.1.0 (Appendix B): 0.1.0 made
  the texts selectable and scrollable (`mdatk/gui/dialogtext.py`), but the cut and the fixed size remain.
- **Where:** `mdatk/gui/app.py:2570`, `:2580` (`format_secondary_text(detail[:3000])` in `error` / `info`),
  `mdatk/gui/rnado.py:336`; `mdatk/gui/app.py:2585-2594` (`show_image`: `screen_w = 1200`, height not
  considered).
- **Description:** a long problem list or traceback loses its end in the dialog without any mark (the log
  has the full text, but nothing says so); a portrait still (e.g. 1080 × 1920) is shown at full height,
  beyond a laptop screen.
- **Expected:** the whole detail in the (already scrollable) box, or a visible "… (full text in the log)";
  the still scaled to fit the monitor's work area in both directions.
- **Suggested fix:** MDBridge shows the whole text, selectable and scrollable, and fits the still to 90 % of
  the available screen (aspect kept). In MDATK: drop `[:3000]` (dialogtext already scrolls long text) and
  scale by `min(1, 0.9 * work_w / w, 0.9 * work_h / h)` from `Gdk.Monitor.get_workarea()`.
- **Test idea:** `win.error("x", "y" * 5000)` with a stubbed `run`; assert the label holds 5 000 characters.

### MDATK-BUG-104 — After *Run cleaning now* the preview keeps showing the raw frame

- **MDBridge ref:** WA-05 · **Area:** GUI · **Severity:** minor · **Status:** Confirmed by reading
- **Affected:** 0.0.4.1 … 0.1.1
- **Where:** `mdatk/gui/app.py:2412-2430` (`run_cleaning`'s `done` shows a dialog only), `:2460-2476`
  (`preview_data` is replaced only by *Load real frame*).
- **Description:** the reference frame loaded before cleaning (raw, with water and ions, unimaged) stays on
  the canvas, although `structure_preview()` would now take the cleaned files.
- **Expected:** the preview reloads the reference frame after cleaning when one was loaded.
- **Suggested fix:** MDBridge loads the reference frame again when cleaning ends (and whenever inputs,
  cleaning or groups change). In MDATK: call `load_preview()` from `run_cleaning`'s `done` when
  `self.preview_data` is set.
- **Test idea:** stub the pipeline; after `run_cleaning`, assert `preview_data["cleaned"]` is true.

### MDATK-BUG-105 — On Windows a save fails when another program holds the file for a moment

- **MDBridge ref:** PRJ-13 · **Area:** project files, parameter logs, cache (Windows) · **Severity:** minor ·
  **Status:** Confirmed by reading (needs Windows)
- **Affected:** 0.0.4.1 … 0.1.1
- **Where:** `mdatk/projectfile.py:353-362` (`_write`: `.tmp` then one `os.replace`), `:519`, `:589`;
  `mdatk/io/cache.py:55`, `:85`; the same pattern in `mdatk/workflows/plan.py:773`.
- **Description:** Windows refuses to replace a file while another program has it (or the new file) open,
  even briefly; a virus scanner, the search indexer or a sync client reading the file just written makes the
  save fail with `PermissionError: [WinError 5] Access is denied` (seen in MDBridge's Windows CI with the
  same pattern). The `.tmp` file is then left next to the project.
- **Expected:** the save waits briefly for the other program; a failed save leaves no `.tmp`.
- **Suggested fix:** MDBridge retries the rename for up to 2 s (10 ms, doubling to 100 ms) while Windows
  reports the file busy (access denied, sharing or lock violation), then fails with the system's message as
  before, and removes `<name>.tmp`. One try on Linux and macOS. In MDATK: a small `replace_retrying(tmp,
  path)` helper used by every atomic write, with `try/finally` removing the `.tmp`.
- **Test idea:** on Windows, hold the target open with a handle sharing read/write/delete for 300 ms and
  assert the save succeeds; hold it without delete sharing past 2 s and assert the error and no `.tmp`.

### MDATK-BUG-106 — One missing Wayland library puts every bundled fallback library ahead of the host's

- **MDBridge ref:** PKG-01 · **Area:** Linux AppImage · **Severity:** minor · **Status:** Confirmed by reading
- **Affected:** 0.0.4.1 … 0.1.1
- **Where:** `packaging/linux/AppRun:16-35` (the loop prepends the whole `_internal/fallback` folder to
  `LD_LIBRARY_PATH` and stops at the first library the host lacks), `packaging/linux/build_appimage.sh:62-75`
  (all `libwayland-*` copies go into that one folder). The comments say the copies are used "only on hosts
  that lack them".
- **Description:** a host with `libwayland-client` but without `libwayland-cursor` (or `-egl`) gets all
  bundled Wayland libraries ahead of its own, mixing the AppImage's build of one with the host's driver
  stack for the others. MDBridge saw the same pattern with a larger fallback set break the host's font and
  GL set-up.
- **Expected:** only the libraries the host lacks come from the AppImage.
- **Suggested fix:** MDBridge puts each fallback copy in a folder of its own
  (`fallback/<name>/<name>`) and adds only the folders of the libraries the host lacks (host check: the
  loader cache, the standard library folders, the caller's `LD_LIBRARY_PATH`); libraries built and shipped
  as one unit (GL's glvnd set) are supplied all or none. For MDATK's three Wayland libraries, per-library
  folders are enough (they come from one package, so all-or-none per that package is also reasonable).
- **Test idea:** a container with `libwayland-client` only; `LD_DEBUG=libs` shows the host's client and
  the bundled cursor / egl libraries.

---------------------------------------------------------------------------------------------------

## 5. First-report items fixed upstream since 2026-11-25

None. The newest `origin/main` (2026-11-26 fetch) is still `1efa0ef` (0.1.1), the commit the first report
was checked against, so all 85 items MDATK-BUG-001 … 085 are present unchanged. (Correction to that
report's Appendix B: AP-18 was only partly fixed in 0.1.0; the remaining part is 103 above.)

MDATK's response to both reports (2026-11-27), including 099 here, is recorded in the first report under
"MDATK's response".

## Appendix: items checked and left out

- **Already in the first report** (not repeated): every AP item of the app-layer review: AP-01 → 001
  (also its `visual.zoom` 0.01 → 0.05 and the rotation sliders' −200 → −180 when another axis moves,
  MDBridge PG-01), AP-02 → 002, AP-03 → 003, AP-04 → 004, AP-05 → 005, AP-06 → 006 (FW-08, picks stored
  absolute, is part of the same design question), AP-07 → 007, AP-08 → 008, AP-10 → 010, AP-13 → 009,
  AP-14 / JB-01 → 068, AP-17 / PRJ-12 → 069; the command-line entries of DIFFERENCES ("No tracebacks" → 028,
  CLI-01 → 065, CLI-02 → 084, CLI-03 → 083, CLI-04 → 030, CLI-08 → 085).
- **Fixed upstream (excluded):** CH-05 (DNA with RNA names, 0.1.0), CH-06 (unknown terminal neighbour,
  0.1.0), CH-07 (four-character SEQRES names, 0.0.17; verified: `LDA5` is written `0DA`), CH-08 (`15G` is
  residue 15, 0.1.0; verified), SEL-19 … SEL-22 and CLI-07 (selection error texts, hints, circular groups, empty
  explanations, 0.1.1), EXP-01 (comma-only export subset, 0.0.10), EXT-01 (two external series of one key,
  0.1.0), TOP-01 (chains from bonds, 0.1.0), NA-01 (L-DNA and CCD names, 0.0.13 – 0.0.17), PRJ-10 (log
  marker inside the settings, 0.1.0), PRJ-11 (control characters in comments, 0.0.13), CMD-05 (cleaned
  files matched by key, 0.1.0), ED-06 / GP-01 (group renames, 0.0.10), WN-11 (pages wider than the side
  bar, 0.1.1 flow rows), WA-04 (*Inspect system* on an empty group: 0.1.1's empty-selection explanation now
  gives `group 'Far' (residues '99000-99010') matches no atoms - the structure's residues are numbered
  1-4095` instead of `index 0 is out of bounds`; verified).
- **Not bugs in MDATK (excluded):** MDBridge-only design choices and conventions: CMD-02, CMD-06 … CMD-08,
  CH-14 (Barnaba), CIT-01 … CIT-06 (citation policy), AN-25 … AN-27 (template and typed Lennard-Jones
  parameters, template charges: accuracy improvements over 0.1.1's documented table; its virtual-site and
  `use_parm7` issues are 021 and 026), FONT-01 … 03, PRJ-07 (a hint), WN-01 … WN-10, WN-12 (Windows path
  separators in texts), WN-13, WN-14, ED-01 … ED-05, ED-07, ED-08, FW-*, PV-*, PG-02 … PG-06, GP-02 …
  GP-05, JB-03, JB-04, WA-01, WA-02, WA-06, WA-08, WA-09, LE-01 … LE-05, LE-07 … LE-10, PKG-02, PKG-03.
  `0.1.1`'s generic Lennard-Jones code (`mdatk/analysis/generic_lj.py`) was read for further bugs; none
  beyond 021 found.
