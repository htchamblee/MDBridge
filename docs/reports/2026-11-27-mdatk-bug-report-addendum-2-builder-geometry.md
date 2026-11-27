# MDATK bug report, addendum 2: structure builder geometry (still present in 0.1.2)

- **Date:** 2026-11-27
- **Follows:** `docs/reports/2026-11-25-mdatk-bug-report.md` (001 … 085) and
  `docs/reports/2026-11-26-mdatk-bug-report-addendum.md` (086 … 106). Ids continue at `MDATK-BUG-107`.
- **Source:** MDBridge's port of the non-Watson-Crick builder (milestone 6, package W2-A) and an independent
  check of its geometry against experimental structures: `docs/reports/2026-11-27-non-wc-geometry.md`
  (the method and the full comparison) and owner decisions D18 and D19 (`docs/specs/m6-plan.md`).
- **MDATK version checked:** 0.1.2, commit `2d93338a125a1ee3b3eeeec8675e7a488afbedb8` (2026-11-27).
  `git diff 1efa0ef 2d93338 -- mdatk/builder.py` changes only `build mutate` (099) and the residue-number
  warning (093), so the builds below are the same in 0.0.4.1, 0.1.1 and 0.1.2. `file:line` is 0.1.2.
- **Status:** *Verified* = measured on MDATK's own builds; *Confirmed by reading* = the cited code does what
  is described.
- **Owner decisions:** the user owns both projects. Items 108, 112, 113 and 114 carry the owner's decisions
  for MDBridge, which MDATK may want to follow; the others are plain defects.

## How the geometry was measured

One level is one base pair or base triple. TWIST and RISE are the screw that maps level i onto level i+1
(rigid fit of the base atoms and C1' of every strand of the level, as `builder.py` defines them); stacking
is the separation of consecutive level base planes. The same measurement was applied to experimental
structures with several consecutive stacked levels of the same pairing (wwPDB entries, CC0). The scripts
are not attached; the method is in the geometry report.

## Summary

9 items: 5 wrong results, 4 minor.

| Id | Severity | Title | Status |
|---|---|---|---|
| MDATK-BUG-107 | wrong results | Hoogsteen, reverse-Hoogsteen, reverse-WC, th and trh helices are unstacked | Verified |
| MDATK-BUG-108 | wrong results | RNA is built in forms RNA is not known to adopt (reverse-WC, reverse-Hoogsteen duplex, purine-motif triplex) | Verified |
| MDATK-BUG-109 | wrong results | The reference fit scores contacts between three levels only | Confirmed by reading |
| MDATK-BUG-110 | wrong results | The non-WC closures are ill-conditioned (about 0.8 Å from the BLAS kernel alone) | Verified |
| MDATK-BUG-111 | minor | Strained links and levels without a base-base H-bond are not reported as such | Verified |
| MDATK-BUG-112 | wrong results | hh levels other than A·A have no H-bond and open the helix | Verified |
| MDATK-BUG-113 | minor | A pyrimidine on the purine strand of a Hoogsteen kind distorts the helix | Verified |
| MDATK-BUG-114 | minor | C+ of C+·G Hoogsteen pairs is never protonated, even with hydrogens | Confirmed by reading |
| MDATK-BUG-115 | minor | A-form builds with hydrogens: H2'(i)-H6/H8(i+1) at 1.57 Å | Verified (not yet compared with experiment) |

## Items

### MDATK-BUG-107 — Hoogsteen, reverse-Hoogsteen, reverse-WC, th and trh helices are unstacked

- **Where:** `mdatk/builder.py:1106` (`_HELIX_START`), `:1576` (`_helix_fit`).
- **Description:** only hh's reference helix is realistic. The others have no stacking: consecutive levels
  sit 4.7-8 Å apart against about 3.3 Å in experiment, with twist and rise far from the experimental values.
  The pairs themselves are sound (H-bonds 2.7-2.9 Å, O3'-P 1.60-1.67 Å).

  | Kind | Experiment: twist / rise, stacking | MDATK |
  |---|---|---|
  | hh | 42.3° / 3.76 Å, 3.50 Å (4JRD) | 41.7° / 3.79 Å, 3.56 Å (realistic) |
  | hoogsteen | 30.5° / 2.91 Å, 3.3 Å (RNA X-ray: 3P22, 4PLX, 6SVS) | 43.6° / 3.36 Å, layers 8.0 Å apart |
  | reverse-hoogsteen | 31° / 3.3 Å, 3.5 Å (DNA only) | 44.5° / 3.34 Å, 4.7 Å |
  | reverse-wc | 36-38° / 3.1-3.3 Å, 3.4 Å (DNA only) | 17.6° / 4.17 Å, 5.7 Å |
  | th (pyrimidine-motif triplex) | 30.4° / 3.01 Å, 3.27 Å (RNA X-ray: 3P22, 4PLX, 6SVS, 7JNH) | 23.7° / 5.41 Å, 5.1 Å |
  | trh (purine-motif triplex) | 30.7° / 3.45 Å, 3.46 Å (DNA only) | 7.3° / 5.84 Å, 6.5 Å |

- **Cause:** the fit has no term that rewards stacking, and its objective has many local minima. From
  `_HELIX_START` the hoogsteen fit stops in an unstacked minimum (cost 170-183; experiment-derived starts
  reach a stacked minimum at cost 25). For th and reverse-hoogsteen the objective itself prefers the
  unstacked helices, so new starting values do not help. `_HELIX_START`'s axis values also depend on the
  nucleotide frame.
- **Suggested fix (what MDBridge does, CH-18):** the hoogsteen duplex and the th triplex take their
  geometry from RNA triplex crystal structures: 60 U·A-U and 5 C+·G-C stacked triples (52 steps) of 3P22
  (KSHV PAN ENE), 4PLX (MALAT1 ENE), 6SVS and 7JNH. Per strand a template nucleotide (the 4JRD sugar
  fitted onto that strand's mean sugar, plus the mean base); the mean triple pose per level type; the mean
  step screw. The reference fit is restrained to the experimental screw, poses, H-bonds and stacking
  within their spreads, with contacts over 12 levels and several starts. MDBridge's uniform 14-level
  helices: hoogsteen 33.3° / 3.02 Å, planes 3.35 Å; th 31.6° / 3.16 Å, 3.34 Å; no close contacts with
  hydrogens. The C++ and its Python reference (`tests/parity/improvements_chem.py`, `_patch_ch18`) are
  under review and will be shared when merged.
- **Test idea:** build 12 or more uniform levels of each kind; assert stacking 3.0-3.5 Å and twist / rise
  within two standard deviations of experiment.

### MDATK-BUG-108 — RNA is built in forms RNA is not known to adopt

- **Where:** `mdatk/builder.py:1046-1058` (`PAIRING_LW`, `TRIPLEX_MOTIFS`); `build_helix` accepts every
  pairing for RNA.
- **Description:** RNA is not known to form parallel reverse-WC duplexes (Pallan et al. 2007; Habibian et
  al. 2018), the antiparallel reverse-Hoogsteen duplex, or the purine-motif triplex (it forms only when all
  three strands are DNA; Semerad & Maher 1994). The experimental geometry for these kinds is DNA. With the
  rigid C3'-endo nucleotide a stacked reverse-WC helix cannot close: a stacking restraint raises the link
  cost from 11 to 150-233.
- **Owner decision (MDBridge D18):** RNA builds of reverse-wc, reverse-hoogsteen and the reverse-hoogsteen
  (purine-motif) triplex are refused, in the library, the command line (error, exit 1, nothing written) and
  the builder window (the items greyed out with the reason). They return for DNA and DNA-RNA hybrids with
  DNA structure building. MDATK already refuses these for DNA (`_dna_only`, `:388`); the same message
  style fits the RNA refusal.

### MDATK-BUG-109 — The reference fit scores contacts between three levels only

- **Where:** `mdatk/builder.py:1531-1535` (`_level_distances`: `levels = np.stack([xyz, screw.apply(xyz,
  1), screw.apply(xyz, 2)])`).
- **Description:** contacts are scored between levels i, i+1 and i+2 only. A reference helix that passes
  can clash from 5 levels on: a refitted reverse-hoogsteen reference that passed this check gave 145 close
  contacts in an AAAAAAAA build (strands interpenetrating, 0.68 Å).
- **Suggested fix (MDBridge):** score contacts over 12 levels in every reference fit (hh included), using
  the pairs within 6 Å refreshed between solves.

### MDATK-BUG-110 — The non-WC closures are ill-conditioned

- **Where:** `mdatk/builder.py:1666` (`_close_levels`).
- **Description:** the strand closure of mixed non-WC and triplex sequences ends in different minima for
  tiny changes in arithmetic: the choice of BLAS kernel alone moves atoms by about 0.8 Å. Results are not
  reproducible across machines.
- **Suggested fix (MDBridge CH-16, owner decision D17):** a continuation closure (the link restraints
  switched on gradually from the uniform helix, converged at each stage) with a smooth link-torsion term;
  for the hoogsteen and th kinds the stacking of each step is also held. Reproducible within 0.001 Å
  across kernels and 1e-8 start perturbations.

### MDATK-BUG-111 — Strained links and levels without a base-base H-bond are not reported as such

- **Where:** `mdatk/builder.py` build notes (link lengths outside 1.5-1.7 Å are the only link geometry
  warned).
- **Description:** links can have C3'-O3'-P angles of 125-127° (hh `UUUUUUUU`: 126.9°; experiment
  119.7 ± 1.1°) without a note. Levels placed as mismatches are warned ("no tHH pair; placed from the A·A
  reference pair"), but the note does not say that the level has no base-base H-bond, and a build whose
  every level is such a mismatch (112, 113) gets no note that the helix is not a realistic structure.
- **Suggested fix (MDBridge):** `WARNING: N base pair(s) without a base-base H-bond (no N/O donor-acceptor
  pair within 3.5 Å): …` for any build with a partner; MDBridge also plans link-angle warnings.

### MDATK-BUG-112 — hh levels other than A·A have no H-bond and open the helix

- **Description:** the hh kind's geometry is the parallel poly(rA) duplex (4JRD) and is realistic only for
  A·A. Measured on 0.1.2: `build_helix("UUUUUUUU", pairing="hh")` has no base-base H-bond in any level
  (closest base N/O across the strands 4.9 Å), C3'-O3'-P up to 126.9°, level planes 4.4-4.9 Å apart in a
  12-mer; in `AAAAUAAAA` the U·U level's closest base N/O is 4.5 Å (2.9 Å in the A·A levels).
- **Owner decision (MDBridge D19):** hh builds accept A·A levels only; other bases are refused.

### MDATK-BUG-113 — A pyrimidine on the purine strand of a Hoogsteen kind distorts the helix

- **Description:** in the hoogsteen duplex and the th triplex, strand 1 carries the Hoogsteen edge, which
  only a purine has. Measured on 0.1.2: `build_helix("UUUUUUUUUUUUUU", pairing="hoogsteen")` builds 14
  levels with no base-base H-bond (closest base N/O 3.7 Å), and in a 12-mer the level planes are 1.6-5.1 Å
  apart. Each level is warned as a mismatch ("a strand-1 pyrimidine has no purine Hoogsteen edge"), and the
  build succeeds.
- **Owner decision (MDBridge D19):** such builds are refused, with a message naming the positions.

### MDATK-BUG-114 — C+ of C+·G Hoogsteen pairs is never protonated

- **Where:** `mdatk/builder.py:574, 621, 973` (notes: "written as plain C (N3 unprotonated)").
- **Description:** with `hydrogens=True` the cytosine of a C+·G Hoogsteen pair or C+·G-C triple stays
  unprotonated; every such level has an acceptor-acceptor N3···N7 contact of 2.79 Å.
- **Owner decision (MDBridge D19):** with hydrogens, these cytosines get the N3 proton and the wwPDB CCD
  code for N3-protonated cytidine, `CH`, with atom names from the CCD definition. Programs that cannot read
  `CH` are dealt with when they are used; the CCD standard comes first.

### MDATK-BUG-115 — A-form builds with hydrogens: H2'(i)-H6/H8(i+1) at 1.57 Å

- **Reproduce:**
  ```python
  from mdatk import builder, pdbtools as pt
  import numpy as np
  t = builder.build_helix("GCAUGCAUGCAU", hydrogens=True)
  r = list(pt.read_pdb(t).residues())
  a = {(x.chain, x.resseq, y.name): np.array(y.xyz) for x in r for y in x.atoms}
  print(np.linalg.norm(a[("A", 3, "H2'")] - a[("A", 4, "H6")]))   # 1.57
  ```
- **Description:** 22 inter-residue H-H pairs under 2.0 Å, all H2'(i)-H6/H8(i+1), 1.57-1.59 Å. The
  heavy-atom model has no contact. This is the Arnott A-RNA geometry with hydrogens placed on it. The
  H2'(i)-H6/H8(i+1) distance is known to be short in A-form RNA; whether 1.57 Å is too short has not yet
  been checked against high-resolution experiment, so this item is listed for checking, not as a
  confirmed defect.
- **Test idea:** compare C2'(i)-C6/C8(i+1) with a sub-1.5 Å RNA duplex crystal structure.
