# Crystal against solution: how far may crystal-derived helix geometry be from solution? (2026-11-28)

The non-Watson-Crick builder's Hoogsteen duplexes and pyrimidine-motif triplexes take their twist and
rise from crystal structures (PDB 3P22, 4PLX, 6SVS, 7JNH; DIFFERENCES CH-18, CH-20). The owner asked
that this geometry be handled "while bearing in mind that its native solution conformation may differ
somewhat from what it adopts under crystallization", and whether that has been quantified.

**No triplex or Hoogsteen duplex has both a crystal and a solution structure of the same sequence in
the PDB.** We searched 381 triplex and Hoogsteen entries and cross-matched them against every X-ray and
NMR nucleic-acid-only entry with chains of at most 40 nucleotides. The near misses were the NMR triplex
1BWG, whose two duplex strands are in crystal 8PJL without the third strand, and the echinomycin-DNA
crystals, whose Hoogsteen pairs the drug induces. The difference is therefore estimated from duplexes
that have been solved both ways.

## Method

Helical twist, rise, inclination and x-displacement of each base-pair step were measured from the
standard base reference frames (Olson et al. 2001, J Mol Biol 313:229). The script was written from the
published definitions and is in `tools/crystal_solution/`, which also lists the files measured and
their SHA-256 sums. Each crystal copy and each NMR model was measured on its own, and then the copies
or models were averaged. "Interior" leaves out the two terminal steps at each end.

## Results (crystal minus NMR; twist in degrees, rise in Å)

| Same molecule | Crystal | NMR | Twist, all steps | Twist, interior | Twist, end steps | Rise, all |
|---|---|---|---|---|---|---|
| purine-rich RNA·DNA hybrid r(GAAGAGAAGC)·d(GCTTCTCTTC) | 404D | 1RRD (10 models, NOE) | −1.5 | −2.6 | −0.1 | +0.16 |
| B-DNA dodecamer CGCGAATTCGCG | 355D, 1BNA | 1NAJ, 1DUF (RDC-refined) | +0.4 | −1.4 | +3.6 | +0.06 |
| A-RNA r(GGCGUGCC) with tandem G·U | 2G3S (5 copies) | 1EKD (1 model, partly modelled) | +1.8 | +2.2 | +1.5 | +0.14 |
| HNA·RNA, *different* sequences and lengths | 2BJ6 (4 copies) | 1EJZ (1 model) | −5.0 | −6.5 | −3.4 | −0.17 |

- **The hybrid is the closest analogue to the builder's purine-rich strands.** In the crystal it is
  more A-like than in solution. The DNA sugars are C3′-endo rather than O4′-endo, the twist is
  1.5-2.6° lower, the x-displacement is 1.4 Å more negative and the inclination is 6° larger. The
  authors' own analysis gives 30.4° against 33.7° (Conn et al. 1999, Table 1). They attribute the A′
  form of the crystal to the crystallisation conditions: high salt and MPD or PEG lower the water
  activity, which favours C3′-endo deoxyriboses. The terminal dC20 is C2′-endo because of packing
  contacts. In an all-RNA triplex every sugar is C3′-endo in both states, so this mechanism should act
  less there (an inference, not a measurement).
- **The direction is not fixed.** The B-DNA dodecamer is about as twisted in the crystal as in
  solution, and the A-RNA is more twisted in the crystal.
- **The 3-5° of Istrate et al. 2019** compares HNA·RNA duplexes of different sequence and length. We
  reproduce it (−5.0°), but there sequence and method cannot be told apart.
- **Ends against centre.** In every crystal the terminal steps are 2-4.5° more twisted than the
  crystal's own centre (404D +3.0°, 2BJ6 +2.0°, 355D / 1BNA +4.5°). The NMR ends match their centres
  (+0.5°, −1.2°, −0.5°). This supports "the ends are distorted by packing".
- **Step-to-step variation is not suppressed in crystals.** Within one structure it is 1.2-2.4 times
  larger in the crystals than in the NMR models. NMR ensembles are smoothed by their restraints and
  force field, so neither method measures thermal fluctuations of the steps directly.

## Conclusion for the builder

The solution conformation of the builder's crystal-derived Hoogsteen and triplex geometry may differ
by about 1-3° of mean twist per step and up to about 0.2 Å of rise, in either direction. That is more
than the builder's own deviation from the crystal mean (uniform helices +1.0° triplex, +1.4° duplex;
CH-18), so the remaining twist offsets are within what crystal and solution differ by. The physical
checks (contacts, H-bonds, link lengths and angles, stacking) stay strict, because they are chemistry
rather than crystal artefacts.

## References (from the RCSB primary-citation records and the papers)

- Conn G.L., Brown T., Leonard G.A. (1999) Nucleic Acids Res 27:555-561, doi 10.1093/nar/27.2.555 (404D).
- Gyi J.I., Lane A.N., Conn G.L., Brown T. (1998) Biochemistry 37:73-80, doi 10.1021/bi9719713 (1RRD, 1AXP).
- Shui X. et al. (1998) Biochemistry 37:8341, doi 10.1021/bi973073c (355D).
- Drew H.R. et al. (1981) PNAS 78:2179, doi 10.1073/pnas.78.4.2179 (1BNA).
- Wu Z., Delaglio F., Tjandra N., Zhurkin V.B., Bax A. (2003) J Biomol NMR 26:297, doi 10.1023/A:1024047103398 (1NAJ).
- Tjandra N., Tate S., Ono A., Kainosho M., Bax A. (2000) J Am Chem Soc 122:6190, doi 10.1021/ja000324n (1DUF).
- Jang S.B. et al. (2006) Biophys J 90:4530, doi 10.1529/biophysj.106.081018 (2G3S).
- Chen X. et al. (2000) Biochemistry 39:8970, doi 10.1021/bi992938e (1EKD).
- Maier T. et al. (2005) J Am Chem Soc 127:2937, doi 10.1021/ja045843v (2BJ6).
- Lescrinier E. et al. (2000) Chem Biol 7:719, doi 10.1016/S1074-5521(00)00017-X (1EJZ).
- Istrate A., Johannsen S., Istrate A., Sigel R.K.O., Leumann C.J. (2019) Nucleic Acids Res 47:4872,
  doi 10.1093/nar/gkz197 (the 3-5° HNA·RNA comparison).
- Olson W.K. et al. (2001) J Mol Biol 313:229 (the standard reference frame).
