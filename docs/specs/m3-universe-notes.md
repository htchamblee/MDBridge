# Universe, atom-group and trajectory behaviours used by the analyses (MDAnalysis 2.11.0)

Clean-room notes for milestone 3 (plan §2(a), third row): which behaviours of MDAnalysis 2.11.0's
Universe, atom groups and trajectory readers that MDATK v0.0.4.1's analyses rely on are already
pinned by the milestone-2 test data, and a behaviour description of the rest. Observed by running
MDAnalysis 2.11.0 with MDBridge's independent tables patched in (`tests/data/README.md`, first
section). Probe data: `tests/data/analysis/universe_cases.json` and the crafted inputs in
`tests/data/analysis/universe_files/`, written by `tests/parity/make_universe_cases.py`.

"Group" below means an ordered list of atom indices that may repeat (an atom group); "the
Universe's bonds" means the list of bonded atom pairs of the whole system.

## 1. Inventory

| Behaviour | Status | Where |
|---|---|---|
| iterating `trajectory[int64 array]`: order, repeats, current frame afterwards | **new** | §2; `iteration` |
| `ts.time` per frame | pinned | `io/trajectories/cases.json` (`time` of every frame, every reader, as `float.hex`), `io/loader/topology_cases.json`; fancy-index iteration yields the same times (§2) |
| `ts.dimensions` None without a cell, all-zero cells | pinned | null `box` / `dimensions` in `io/trajectories/cases.json` (no-cell NetCDF, mdcrd, rst7, DCD, XTC, TRR; `dcd_zero_cell`), topology dumps (`cryst_zero.pdb`, `zero_box.gro`); rule in `m2-triclinic-box.md` §3 |
| cells with **some** zero lengths or angles | **new** | §3; `dimensions` |
| whether the atoms have charges (`hasattr`) per format | presence pinned, error **new** | dumps: `charges` null / non-null; §4 |
| masses per format | pinned | dumps `masses` (parm7: the file's MASS section; PDB / GRO: guessed from the type with the independent tables and the residue rule); always f64, never absent |
| the Universe's bonds of formats without bonds (error type, message) | **new** | §5; dumps only give `bonds: null` |
| a group's bonds: membership, order; the order of the Universe's bonds | **new** | §5 (dumps sort the pairs themselves) |
| positions dtype | pinned | f32 everywhere (`positions` of the dumps and trajectory cases); recorded again per topology |
| a group's residues; the atoms of a residue group | **new** | §6 |
| fragment indices | **new** | §7 |
| chain IDs: presence per format | pinned | dumps `chainIDs` (PDB: always; parm7: only with a used `RESIDUE_CHAINID` section, `io/parm7/topology/chain_*`; GRO: never); error **new**, §4 |
| boolean-mask subsets, union of groups | **new** (small) | §6 |

## 2. Iterating over frames given by an index array

MDATK's frame loop iterates the trajectory indexed by a one-dimensional int64 array of frame
numbers (its validated frame list), after every analysis' setup; the rmsd analysis' setup first
moves the trajectory to the reference frame to read its positions.

Observed for NetCDF, XTC, DCD, TRR, multi-model PDB and GRO readers alike:

* **Indexing** checks every index at once, before any frame is read. A negative index k counts
  from the end (k + n_frames). An index ≥ n_frames, or < −n_frames, raises `IndexError` with the
  message `Index K exceeds length of trajectory (N).`, where K is the index **after** adding
  n_frames to a negative one (index −4 of a 3-frame trajectory: `Index -1 exceeds length of
  trajectory (3).`) and N is n_frames. Nothing is read and the current frame is unchanged.
* **Iteration** yields the frames in array order, repeats included (`[2, 2, 0, 2, 0, 0]` yields
  six frames), each exactly what reading that frame directly gives (positions, cell, time).
* **After the iteration has run to its end** the trajectory is back at **frame 0** (frame 0 read
  again: its positions, cell and time are current), whatever the last yielded frame was and
  whatever frame was current before. This also holds for an **empty** array: nothing is yielded
  and the current frame becomes 0.
* **If the consumer stops early** (breaks out of the loop), the current frame is the last frame
  yielded (no rewind).
* A seek before the loop has no influence on what is yielded.

Consequence for MDATK: after `run_analyses` the trajectory is at frame 0. During setup, the
current frame is the one the previous setups left: frame 0 after loading, but **the rmsd reference
frame for every analysis set up after rmsd** (analyses are set up in turn before the loop). The
"atoms" view's setup guesses bonds from the current positions when the topology has none
(`atom_bonds`), so when it runs after rmsd it uses the rmsd reference frame, not the first frame.
A port must keep the same "current frame" state across setups.

(Not used by MDATK: an int32 index array fails on the NetCDF reader with `TypeError: frame must be
a positive integer or zero`; boolean arrays select frames.)

## 3. Cells with zero lengths or angles

The "no cell" rule (a frame has no cell when its three lengths are all exactly zero) is pinned. The
crafted files show what the readers do with cells that are only partly zero, written by
MDAnalysis' own writers from the cells below:

| Cell written | NetCDF, DCD | XTC, TRR |
|---|---|---|
| (0, 20, 20, 90, 90, 90) | (0, 20, 20, 90, 90, 90) | none |
| (0, 0, 0, 90, 90, 90) | none | none |
| (0, 0, 0, 0, 0, 0) | none | none |
| (20, 20, 20, 0, 90, 90) | (20, 20, 20, 0, 90, 90) | none |
| (20, 21, 22, 90, 90, 90) | as written | (20, 21.000002, 22, 90, 90, 90) |

NetCDF and DCD store lengths and angles and return a partly zero cell as it is (it is a cell, with
zeros in it); XTC and TRR store box vectors, and the conversion of `m2-triclinic-box.md` makes such
a cell invalid, so it reads as no cell. MDATK's own check (positive lengths) then decides whether
the cell is usable; that check is MDATK code.

## 4. Optional attributes: charges, chain IDs

MDATK asks `hasattr(atoms, "charges")` (energy) and wraps chain-ID access in a handler that catches
every exception (structure). When an attribute is absent, reading it raises `NoDataError`, an
exception that **is both an AttributeError and a ValueError** (so `hasattr` is false, and the
handlers catch it). Messages:

* charges: `This Universe does not contain charge information`
* chain IDs: `This Universe does not contain chainID information`
* positions of a Universe without coordinates (a parm7 file alone):
  `This Universe has no coordinates`

Per format (the dumps pin the same presence): parm7 — charges always (f64), chain IDs only with a
used `RESIDUE_CHAINID` section; PDB — never charges (formal charges are a different attribute and
do not count), chain IDs always (values pinned by the dumps); GRO — neither.

## 5. Bonds

**Formats without bonds.** GRO files, PDB files without CONECT records, and PDB files whose CONECT
records name no partner (`universe_files/conect_lonely.pdb`: `CONECT    1` alone) have **no bond
attribute at all**, not an empty one. Reading the Universe's bonds, a group's bonds, its fragment
indices or its fragments then raises `NoDataError` (AttributeError and ValueError, §4):

* bonds (of the Universe or of any group, the empty group included):
  `This Universe does not contain bonds information`
* fragment indices: `AtomGroup.fragindices not available; this requires bonds`
* fragments: `AtomGroup.fragments not available; this requires bonds`

Every parm7 file probed has bonds (an empty bond section still leaves the other section's bonds); a PDB file with at least one valid CONECT pair has bonds (which
pairs is pinned by the loader dumps).

**Order of the Universe's bonds.** Each bond is an index pair (i, j) with i < j, each pair once;
the list is sorted by i, then j — whatever the order of the CONECT records or parm7 bond sections
(`frag_conect.pdb` lists them out of order and in both directions; demo: 828 pairs, sorted).

**A group's bonds** are the bonds of the Universe with **at least one** atom in the group (the
other atom may be outside it), each bond once however often its atoms occur in the group, in the
Universe's order (sorted as above), **not** in the group's order. Examples (demo):

| Group | Its bonds |
|---|---|
| [30, 5, 5, 2] | (1,2) (2,3) (4,5) (5,6) (5,11) (19,30) |
| [400, 30, 31] | (19,30) (20,31) (399,400) (400,404) |
| [772, 0, 771] (two Mg ions and atom 0) | (0,1) (0,8) (0,15) |
| [771] (an Mg ion) | none (empty list, no error) |
| [] | none |

A group of a Universe with bonds but without bonds of its own gives an empty list (length 0), not
an error; MDATK's hbonds treats "length > 0" as "bonds known".

## 6. Residues, masks, unions

* A group's **residues**: the distinct residues of its atoms, in increasing residue index, each
  once (group [431, 30, 31] → residues 0, 13; [772, 0, 771] → 0, 25, 26).
* The **atoms of a residue group**: for each residue of the residue group in **its** order,
  repeats included, that residue's atoms in increasing index. Residue group [3, 1, 1] gives the
  atoms of residue 3, then those of residue 1 twice; [26, 0] gives atom 772, then 0–33. (For the
  residues of an atom group, which are sorted and unique, this is all their atoms in index order.)
* A **boolean-mask** subset keeps the selected entries in the group's order, repeats included
  ([30, 5, 5, 2] with mask [T, T, T, F] → [30, 5, 5]).
* The **union** `g | h` is the sorted set of the atoms of both, without repeats ([5, 3] | [3, 1, 1]
  → [1, 3, 5]).

## 7. Fragment indices

A fragment is a connected component of the bond graph over all atoms (an atom without bonds is a
fragment of its own). Fragments are numbered 0, 1, 2, … in increasing order of their **smallest
atom index**. A group's fragment indices are the Universe-wide indices of its atoms, in the
group's order, repeats included. A group's fragments are the distinct fragments of its atoms in
increasing fragment index, each as all of its atoms (also those outside the group) in index order.

Examples: `frag_conect.pdb` (8 atoms; bonds 1–2, 1–6, 3–4, 4–7) → fragment indices
[0, 1, 1, 2, 2, 3, 1, 2]; group [7, 4, 1] → [2, 2, 1], fragments [1, 2, 6], [3, 4, 7].
Demo: five fragments (the two strands, then each Mg ion). `io/parm7/topology/chain_each.parm7`
(bond 0–39): atoms 0, 1, 39 are fragment 0, and the later numbers follow the smallest atom of each
fragment ([0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 6, 7, 5, 8, …, 0]). Without bonds: the error of §5.

## 8. The probe data (`tests/data/analysis/universe_cases.json`)

* `iteration`: per system (`topology`, `trajectory`, paths relative to `tests/data`; the demo
  NetCDF / XTC / DCD / TRR, `io/loader/files/models_npt.pdb`, `io/loader/files/seven.gro`) and
  scenario: `seek` (frame moved to first, or null), `indices`, `stop_after` (frames consumed
  before stopping; null = all), `yielded` (`[frame, time]` in order), `after` (current frame and
  time afterwards), `error` / `at` (`"index"`: raised when indexing, `"iterate"`: while iterating;
  null otherwise). Times are `float.hex` of the frame time.
* `dimensions`: per crafted file the cell written and the cell read for each frame (null = none).
* `topology`: per file (`analysis/universe_files/...` is relative to `tests/data`):
  `has_charges`, `has_chainIDs` (hasattr), `charges` / `chainIDs` / `positions_dtype` /
  `bonds` / `bonds_len` / `fragindices` / `n_fragments` as `{value, error}`, `masses_dtype`, and
  `groups`: per index list its `bonds`, `fragindices`, `fragments` (`{value, error}`), `residues`,
  `residue_atoms`, a `mask` with the `masked` result, and `union_with` / `union`.
* `residues`: residue groups given as index lists and their atoms (`atoms`, `atom_resindices`).
* Errors are `{type, message, attribute_error, value_error}`.
