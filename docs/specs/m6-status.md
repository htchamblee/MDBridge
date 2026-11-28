# Milestone 6: status (2026-11-28)

Milestone 6 is complete except the other pairings and triplexes of the structure builder (W2-A),
which wait for their integration <!-- W2-A -->: every other package is done, and CI is green on the
branch head before W3-C (d15a295). W3-C changes documents only.

Plan: `docs/specs/m6-plan.md`; owner decisions D1-D14 accepted (2026-11-26), D15 (CHD-01: guanine
amino hydrogens follow the wwPDB CCD), D16 (CHD-02: non-WC builds in the fitted 4JRD frame), D17
(CH-16: a continuation closure for the non-WC builds, 2026-11-27) and D18 (CH-18: realistic non-WC
geometry; RNA builds of the kinds RNA does not form refused, 2026-11-27). <!-- W2-A: D19 -->

Reports kept for reference (`docs/reports/`):

- `2026-11-26-mdatk-bug-report-addendum.md`: 21 new items still present in MDATK 0.1.1, among them
  the chemistry items 087-091, 093-096 and 099 (DIFFERENCES CH-01 ... CH-04, CH-09 ... CH-13 and
  CH-15).
- MDATK's response to both reports (0.1.2, 2026-11-27), recorded in `2026-11-25-mdatk-bug-report.md`:
  the reports' chemistry items are fixed in 0.1.2, 099 after the owner's decision D11. MDBridge adopted 0.1.2's
  fixes of items 008, 018, 029, 069, 071 and 077 (1b032c5); the DIFFERENCES entries say how 0.1.2
  compares.
- `2026-11-27-non-wc-geometry.md`: MDATK's non-WC reference helices measured against experimental
  structures (D18).
- `2026-11-27-mdatk-bug-report-addendum-2-builder-geometry.md`: the builder geometry items 107-115,
  checked on MDATK 0.1.2.

## Wave 0: done (CI green at a7cacc5)

Chem generators and `improvements_chem.py` (128 command-line cases pending until the commands
existed), data sources (CCD 0A/0C/0G/0U, PDB 4JRD, Arnott A-RNA; `gen_chem_tables.py`), numerics
library (`src/numerics`: ndimage and densityCentroid moved from render, vendored ducc0 FFT,
Powell/Brent, least_squares TRF with LSMR, pair search; NU-01/02). Generators pin numpy/OpenBLAS
kernels; numerics tests are libm-aware (bit-exact with glibc).

Known fragile result for W2-A, decided since by D17 and D18: MDATK's non-WC / triplex least-squares
closures are ill-conditioned (~0.8 Å spread from the BLAS kernel choice alone, far beyond D10's
tolerance).

## Wave 1: done (CI green at 8d75a4a)

PDB hygiene (`chem/pdbtools.hpp`), annotation (`chem/annotate.hpp`), dictionaries and SMILES
(`chem/dictionaries.hpp`, `chem/smiles.hpp`), the Watson-Crick helix builder (`chem/helix.hpp`) and
the density fit (`chem/density_fit.hpp`); each matches its recorded MDATK cases (byte for byte or
within D10). The regenerate job has a kernel probe: on runners whose pinned kernels give other bits,
MDATK's path-dependent chem results (density-fit evaluation counts, non-WC closures) are compared
structurally (interim; tighten once W2-A makes the non-WC builds well-conditioned).

## Wave 2: W2-B and W2-C done (CI green at 89a7f0e); W2-A pending integration <!-- W2-A -->

W2-B mutations and fit-to-map, W2-C the `pdb`, `build`, `dict` and `density fit` / `simulate`
commands, the pipeline's automatic fit and `chem/generator_api.hpp` (merged at bf1f67d; two test bugs
the sanitizer job found fixed at c4342d6).

<!-- W2-A: replace this paragraph with W2-A's summary when it is integrated. -->
W2-A (other pairings and triplexes) is pending integration: it implements D17 (CH-16, MDATK's
`_close_levels` replaced by a continuation closure) and D18 (CH-18) with their expected data. Until
then the library answers the other pairings with `pairing ... is built by the fitted helix builder
(not available yet)` (TW-07), and the 5 pairing cases of the command line are pending
(`cli_parity.py` PENDING).

## Wave 3: W3-A, W3-B and W3-C done (W3-A / W3-B: CI green at 5f404b0)

- W3-A, the density fit in the window: *Global search* and *Refine* on the Structure page run the
  command line's fit as a window job, the placement stored as one undo step, nothing on Cancel, a map
  changed during the fit not given the old placement, Render / Render still keep an automatic fit's
  placement (DIFFERENCES WF-01 ... WF-03; GUI case `canvas/density_fit`).
- W3-B, the PDB clean-up and structure builder windows behind the host interface (D1), opened from
  the Workflow menu (DIFFERENCES TW-01 ... TW-08; GUI cases `tools/*`; ctest `gui.tools`,
  `gui.tools.scale150`, `gui.tools.scale200`). The non-WC builds of the builder window wait for W2-A
  (`ToolTests::builderNonWcHelixBuildsInTheBackground` is marked `PENDING(W2-A)`). <!-- W2-A -->
- Integrated together: a GUI case's `fitted` (`{file: tolerance}`, need `fitted`) is one mechanism
  for both (tests/data/README.md "gui"): fitted files line by line with their decimals within the last
  printed digit or the tolerance, a project's `[density]` table within D10's angle / shift / CC
  tolerances, and the case's log, snapshots and reports with their decimals within the last digit.
  Two test bugs CI found fixed at 5f404b0.
- W3-C, the docs, reviewed against the code: README (what milestone 6 adds, checked against the
  program's `--help`), BUILDING (the commands without Qt, the tool-window tests, the chemistry
  steps of the regenerate job, the kernel pinning, the data-source extractors, ducc0's check),
  LICENSING (no Python, NumPy, SciPy or MDAnalysis bundled; the NumPy summation in the density fit),
  THIRD_PARTY_NOTICES (the same, and the chemistry data), ARCHITECTURE (status, libraries, module
  dependencies, the tools host, the generator API) and DIFFERENCES "Chemistry" (every CH item with
  what, why, the MDATK version compared and its C++ and command-line tests; the citation keys of the
  chem commands; DIC-01 / DIC-02 checked on 0.1.2).

Pending for the milestone: W2-A's integration, and with it its DIFFERENCES entries (CH-16 ...),
TW-07, the README's and ARCHITECTURE's status lines and the parts of this file marked W2-A.
<!-- W2-A -->
