# Non-Watson-Crick helix geometry against experimental structures (2026-11-27)

This report checks the reference helices of the structure builder's non-Watson-Crick kinds against
experimental structures. It covers MDATK 0.0.4.1, whose builder code is the same in 0.1.1, and
MDBridge's W2-A port with the CH-17 refit.

How the comparison was made:

- One level is one base pair (or one base triple). TWIST and RISE are the screw that maps level i onto
  level i+1: a rigid fit of the base atoms and C1' of every strand in the level, as `builder.py`
  defines them.
- The same measurement was applied to experimental structures that have several consecutive stacked
  levels of the same pairing, and to the uniform helices rebuilt from each version's fitted reference.
- Each family was then re-measured by a second, independent check that tried to refute the result.
  The reverse-WC result had no dedicated second check. An independent check of another family
  measured it in passing and also found it unstacked.

| Kind | Experiment: twist / rise, stacking | MDATK | CH-17 refit | Verdict |
|---|---|---|---|---|
| hh | 42.3° / 3.76 Å, 3.50 Å (4JRD; also 5VXQ, 5K8H, 5TGG) | 41.7° / 3.79 Å, 3.56 Å | same | realistic |
| hoogsteen | 30.5° / 2.91 Å, 3.3 Å (RNA X-ray: 3P22, 4PLX, 6SVS; 55 steps) | 43.6° / 3.36 Å, layers 8.0 Å apart | 41.2° / 3.55 Å, 7.9 Å | unstacked in both |
| reverse-hoogsteen | 31° / 3.3 Å, 3.5 Å (DNA only: 136D, 8PWR NMR; 431D, 272D, 3L1Q X-ray) | 44.5° / 3.34 Å, 4.7 Å | 69.5° / 3.44 Å, 7.3 Å; strands interpenetrate from 5 levels (0.68 Å) | unstacked in both; CH-17 refit impossible |
| reverse-wc | 36–38° / 3.1–3.3 Å, 3.4 Å (DNA only: 1JUU, 1BE5, 1R2L) | 17.6° / 4.17 Å, 5.7 Å | same | unstacked in both |
| th (pyrimidine-motif triplex) | 30.4° / 3.01 Å per triple, 3.27 Å (RNA X-ray: 3P22, 4PLX, 6SVS, 7JNH) | 23.7° / 5.41 Å, 5.1 Å | 23.7° / 5.43 Å, 5.1 Å | unstacked in both |
| trh (purine-motif triplex) | 30.7° / 3.45 Å, 3.46 Å (DNA only: 134D, 135D, 136D, 8PWR) | 7.3° / 5.84 Å, 6.5 Å | 16.8° / 5.80 Å, 5.5 Å | unstacked in both |

In every kind the pairs themselves are sound: Hoogsteen and Watson-Crick H-bonds of 2.7-2.9 Å,
O3'-P links of 1.60-1.67 Å.

## Causes

- **No stacking term, many local minima.** The reference fit has no term that rewards stacking, and
  its objective has many local minima.
  - Started from `_HELIX_START`, the hoogsteen fit stops in an unstacked minimum. Experiment-derived
    starts reach a much lower, stacked minimum: cost 25.3 against 170-183, at 39.1° / 2.99 Å.
  - For th and reverse-hoogsteen, the objective itself prefers the unstacked helices. Changing the
    starting values does not fix them.
- **The nucleotide.** The rigid C3'-endo 4JRD nucleotide cannot form a stacked parallel reverse-WC
  helix: with a stacking restraint, the link cost rises from 11 to 150-233. This matches the
  literature: RNA is not known to form parallel reverse-WC duplexes (Pallan et al. 2007; Habibian et
  al. 2018). Purine-motif triplexes form only when all three strands are DNA (Semerad & Maher 1994).
- **Contact check too short.** The fit scores contacts between levels 0-2 only, in MDATK and in the
  port. The CH-17 reverse-hoogsteen helix passes it, although it clashes from 5 levels on: MDATK
  reports 145 close contacts for AAAAAAAA built with that reference.
- **Frame dependence.** `_HELIX_START`'s axis parameters depend on the frame. The CHD-02 frame (D16)
  therefore moved some fitted references away from what MDATK 0.0.4.1 prints, for example
  reverse-hoogsteen from 48.1° / 4.39 Å to 44.5° / 3.34 Å. The new and old values are unstacked alike.
- **Non-A hh pairs.** In hh builds of sequences other than A, the U·U and C·C levels have no base-base
  H-bond: the closest polar atoms are 4.9-5.3 Å apart.

## Owner decision

The decision is recorded as D18 in `docs/specs/m6-plan.md`.

- The geometry is fixed before W2-A merges.
- The hoogsteen duplex and the th triplex take their geometry from RNA triplex structures: viral
  and cellular elements such as the KSHV PAN ENE (3P22), the MALAT1 ENE (4PLX), 6SVS, 7JNH and the
  telomerase pseudoknots.
- RNA builds of reverse-WC, reverse-hoogsteen and trh are refused for now. They will be built for DNA
  and DNA-RNA hybrids together with DNA structure building.

The scripts and measurements are in the session scratchpad (`geomcheck/`) and are not part of the
repository. The coordinates come from wwPDB entries (CC0).

## Final geometry (owner decision D19, 2026-11-27, evening)

A review of the CH-18 builds found the uniform A and G helices sound but the mixed sequences strained
(C3'-O3'-P up to 127.9°, stacking 0.2 Å loose, the Hoogsteen duplex's steps alternating 28.3° / 38.8°),
the reference fits 2.7° (duplex) and 0.12 Å (triplex) off their screw, a 3.00 Å contact between
stacked bases and statistics that counted copies of one molecule as independent. D19 (docs/specs/
m6-plan.md) and docs/DIFFERENCES.md CH-18 (revised), CH-19 ... CH-21 record what changed:

- the strand-1 pyrimidines of the Hoogsteen kinds and hh's non-A levels are refused (CH-19); the C+
  of the G·C+ Hoogsteen pairs is the CCD's N3-protonated cytidine CH (CH-21);
- the statistics count each molecule once (24 independent steps of 53; CH-20), and the selection of
  the triples is a committed tool (`tools/select_pdb_triplexes.py`; it added 4PLX C75·C38-C15);
- the reference fit is restrained to the experimental means by their standard errors, the contact
  walls follow the report's limits (0.1 Å outside), stacked bases are kept from H-bonding, and the
  closure moves every partner within its level with its pairs' H-bonds, the steps and the levels'
  tilt restrained to the experiment.

**The numbers of this report's first section.** The table at the top (30.5° / 2.91 Å over 55 steps for
the Hoogsteen pairs, 30.4° / 3.01 Å for the triples) was measured on the check's own selection, which
counted a pair as Hoogsteen-paired with a looser criterion (it included the 6SVS steps A76-A78 and
A78-A107, whose Hoogsteen N6-O4 is 3.76 Å, or whose residue is modified) and every copy of a molecule.
The builder's numbers are those of `data/sources/pdb/triplex_levels.tsv` (53 steps; the means over
the 24 independent steps: duplex 30.48° / 3.007 Å, triplex 30.40° / 3.021 Å, base planes 3.270 / 3.263
Å; CH-18's table had 30.62° / 2.996 Å and 30.52° / 3.043 Å over all 52 steps of its selection).

**Crystal and solution.** The geometry is crystal-derived (X-ray, 2.5-3.1 Å). The NMR structures
differ: the telomerase pseudoknot 2M8K 33.0° / 3.36 Å with base planes 3.48 Å apart, the hTR
pseudoknots 28.7° / 3.72 Å. The solution conformation may differ by a few degrees of twist per step (to
be quantified). Leaving out the end-of-run steps or the steps
with a lattice contact moves no mean by its spread (data/sources/pdb/README.md), so every step is used.
The physical checks (contacts, H-bonds, link lengths and angles, stacking) are held strictly; twist
and rise are crystal values.

**Final numbers** (C++ builds with hydrogens, the measures of this report; every step of every model;
experiment: the steps of the table, minimum-maximum, the means over the independent steps):

{{FINALTABLE}}

Still outside the experiment: the purine strand's stacking reaches 3.48 Å in the mixed Hoogsteen
duplexes (experimental maximum 3.42 Å) and the Watson-Crick strand's 3.81 Å in the mixed triplex
(3.79 Å); the inclination of one triplex step reaches 14.7° (the triples' maximum 14.2°, the
Hoogsteen pairs' 17.1°); the pyrimidine strands stack 3.5-3.6 Å on average against 3.41-3.42 Å. The
uniform helices' twist stays 1.0° (triplex) and 1.4° (duplex) above the experimental mean: without the
backbone torsion terms the restrained fit sits on the mean screw (30.6-30.7°), but the rigid C3'-endo
4JRD backbone with MDATK's staggered torsion preferences then needs eclipsed torsions or clashes; the
offset is within the steps' spread (1.85° / 3.06°) and within the expected crystal-solution difference.

The scripts and measurements of this check are in the session scratchpad (`w2af/meas.py` and its
outputs), not in the repository.
