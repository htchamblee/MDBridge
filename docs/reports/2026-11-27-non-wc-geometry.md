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

| model | twist ° | rise Å | base planes Å | inclination ° | strand stacking Å | C3'-O3'-P ° | O3'-P-O5' ° | epsilon ° | pair N/O Å | stacked N/O min Å | contacts |
|---|---|---|---|---|---|---|---|---|---|---|---|
| experiment, Hoogsteen pairs (53 steps) | 30.7 (21.6-37.7) | 2.96 (2.55-3.93) | 3.28 (3.06-3.49) | 10.3 (3.0-17.1) | 3.21 (≤3.42) / 3.42 (≤3.71) | 119.8 (117.1-124.6) | 103.1 (100.2-106.4) | -188 … -101 | 2.80 (2.57-3.02) | 2.82 | - |
| experiment, triples (53 steps) | 30.6 (26.9-34.9) | 3.04 (2.78-3.78) | 3.27 (3.13-3.49) | 9.0 (0.7-14.2) | 3.21 (≤3.42) / 3.41 (≤3.79) / 3.42 (≤3.71) | 119.7 (117.1-124.6) | 103.1 (100.2-106.4) | -188 … -101 | 2.81 (2.41-3.17) | 2.82 | - |
| hoA14 | 31.9 (31.9-31.9) | 3.02 (3.02-3.02) | 3.29 (3.29-3.29) | 6.9 (6.9-6.9) | 3.29 (≤3.29) / 3.57 (≤3.57) | 119.9 (117.9-121.9) | 103.7 (102.3-105.1) | -161 … -149 | 2.83 (2.83-2.83) | 3.28 | 0 |
| hoG14 | 32.3 (27.6-34.2) | 2.94 (2.79-3.07) | 3.31 (3.28-3.34) | 11.4 (9.3-14.5) | 3.15 (≤3.24) / 3.52 (≤3.58) | 120.0 (119.3-120.5) | 104.2 (103.9-104.4) | -163 … -152 | 2.83 (2.79-2.98) | 3.10 | 0 |
| hoGA14 | 32.6 (30.1-35.5) | 2.97 (2.77-3.08) | 3.30 (3.28-3.32) | 9.6 (7.3-12.8) | 3.20 (≤3.40) / 3.54 (≤3.59) | 120.1 (119.6-120.6) | 104.2 (103.9-104.5) | -164 … -156 | 2.83 (2.78-2.97) | 3.15 | 0 |
| hoMix10 | 32.3 (29.2-33.7) | 2.97 (2.84-3.13) | 3.30 (3.28-3.34) | 10.2 (7.3-13.0) | 3.20 (≤3.47) / 3.55 (≤3.64) | 119.9 (119.4-120.4) | 104.2 (104.0-104.4) | -164 … -154 | 2.83 (2.79-2.86) | 3.10 | 0 |
| hoMix32 | 32.1 (28.2-35.0) | 2.98 (2.74-3.19) | 3.32 (3.28-3.35) | 10.1 (7.2-14.0) | 3.22 (≤3.42) / 3.58 (≤3.75) | 120.0 (119.4-120.6) | 104.2 (103.9-104.5) | -167 … -150 | 2.85 (2.79-2.99) | 3.10 | 0 |
| thA14 | 31.4 (31.4-31.4) | 3.05 (3.05-3.05) | 3.28 (3.28-3.28) | 7.5 (7.5-7.5) | 3.31 (≤3.31) / 3.54 (≤3.54) / 3.58 (≤3.58) | 119.5 (117.1-121.4) | 103.6 (101.7-104.7) | -161 … -148 | 2.84 (2.83-2.86) | 3.14 | 0 |
| thG14 | 32.0 (29.1-35.6) | 3.04 (2.82-3.19) | 3.30 (3.28-3.31) | 11.1 (8.6-17.6) | 3.23 (≤3.40) / 3.56 (≤3.71) / 3.58 (≤3.72) | 119.8 (119.6-120.6) | 104.1 (103.9-104.6) | -164 … -153 | 2.83 (2.79-3.05) | 3.10 | 0 |
| thMix10 | 31.6 (29.8-32.8) | 3.05 (2.83-3.17) | 3.30 (3.28-3.31) | 9.1 (5.7-14.7) | 3.26 (≤3.43) / 3.59 (≤3.81) / 3.60 (≤3.69) | 119.9 (119.6-120.6) | 104.2 (103.9-104.5) | -163 … -154 | 2.84 (2.79-3.01) | 3.10 | 0 |
| thMix32 | 31.5 (28.3-34.5) | 3.04 (2.73-3.24) | 3.30 (3.28-3.32) | 10.8 (6.1-17.5) | 3.27 (≤3.46) / 3.60 (≤3.85) / 3.61 (≤3.79) | 119.9 (119.6-120.7) | 104.2 (103.9-104.5) | -171 … -149 | 2.84 (2.79-3.03) | 3.10 | 0 |
| hhA14 | 41.7 (41.7-41.8) | 3.79 (3.79-3.79) | 3.79 (3.79-3.79) | 0.0 (0.0-0.0) | 3.59 (≤3.59) / 3.59 (≤3.59) | 120.7 (120.6-120.7) | 103.4 (103.4-103.4) | -146 … -145 | 2.74 (2.74-2.74) | 3.66 | 0 |

(ho: the Hoogsteen duplex, th: the pyrimidine-motif triplex; A14 / G14: 14 levels of A / G; GA14:
GAGA...; Mix10: AAGAGAAGGA; Mix32: AGGAAGAGAAGGAAGAGGAAAGAGAAGGAGAA. Means with the range over the
steps (or links, pairs); strand stacking per strand, mean and maximum; "pair N/O": the closest
donor-acceptor pair of a level's paired bases; "stacked N/O": the closest N/O pair of bases at
consecutive levels; contacts: the builder report's limits, hydrogens included.)

Still outside the experiment: the purine strand's stacking reaches 3.47 Å in the mixed Hoogsteen
duplex (experimental maximum 3.42 Å), the triplex's Watson-Crick and Hoogsteen strands 3.85 and 3.79 Å
in the 32-level mixed triplex (3.79 / 3.71 Å); the inclination of a triplex step reaches 14.7° in
AAGAGAAGGA, 17.5° in the 32-level mixed triplex and 17.6° in G14 (the triples' maximum 14.2°, the
Hoogsteen pairs' 17.1°); the triplex G14 has one step of 35.6° (34.9°) and the 32-level mixed triplex
one rise of 2.73 Å (2.78 Å); the pyrimidine strands stack 3.5-3.6 Å on average against 3.41-3.42 Å.
These are single steps at C+·G-C levels and helix ends; the cause is not isolated (the tilt
restraint's spread, 1.9°, does not keep them; a tighter one was not tried). The
uniform helices' twist stays 1.0° (triplex) and 1.4° (duplex) above the experimental mean: without the
backbone torsion terms the restrained fit sits on the mean screw (30.6-30.7°), but the rigid C3'-endo
4JRD backbone with MDATK's staggered torsion preferences then needs eclipsed torsions or clashes; the
offset is within the steps' spread (1.85° / 3.06°) and within the expected crystal-solution difference.

The scripts and measurements of this check are in the session scratchpad (`w2af/meas.py` and its
outputs), not in the repository.
