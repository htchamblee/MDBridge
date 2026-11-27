# Milestone 6: status

Plan: `docs/specs/m6-plan.md`; owner decisions D1-D14 accepted (2026-11-26), D15 (CHD-01: guanine
amino hydrogens follow the wwPDB CCD), D16 (CHD-02: non-WC builds in the fitted 4JRD frame) and D17
(CH-16: a continuation closure for the non-WC builds, 2026-11-27).

## Wave 0: done (CI green at a7cacc5)

Chem generators and `improvements_chem.py` (128 command-line cases pending until the commands
exist), data sources (CCD 0A/0C/0G/0U, PDB 4JRD, Arnott A-RNA; `gen_chem_tables.py`), numerics
library (`src/numerics`: ndimage and densityCentroid moved from render, vendored ducc0 FFT,
Powell/Brent, least_squares TRF with LSMR, pair search; NU-01/02). Generators pin numpy/OpenBLAS
kernels; numerics tests are libm-aware (bit-exact with glibc).

Known fragile result for W2-A: MDATK's non-WC / triplex least-squares closures are ill-conditioned
(~0.8 Å spread from the BLAS kernel choice alone, far beyond D10's tolerance). Treat it as a
fragile result to fix, not to copy; options go to the owner before W2-A finishes.

## Wave 1: done (CI green at 8d75a4a)

PDB hygiene (`chem/pdbtools.hpp`), annotation (`chem/annotate.hpp`), dictionaries and SMILES
(`chem/dictionaries.hpp`, `chem/smiles.hpp`), the Watson-Crick helix builder (`chem/helix.hpp`) and
the density fit (`chem/density_fit.hpp`); each matches its recorded MDATK cases (byte for byte or
within D10). The regenerate job has a kernel probe: on runners whose pinned kernels give other bits,
MDATK's path-dependent chem results (density-fit evaluation counts, non-WC closures) are compared
structurally (interim; tighten once W2-A makes the non-WC builds well-conditioned).

## Wave 2: W2-B and W2-C done; W2-A pending

W2-B mutations and fit-to-map, W2-C the `pdb`, `build` and `density` commands and the pipeline's
automatic fit (merged at bf1f67d). W2-A (other pairings and triplexes) is pending: the owner decided
the conditioning question with D17 (CH-16, MDATK's `_close_levels` replaced by a continuation
closure); the package implements it and its expected data. Until then the library answers the
other pairings with `pairing ... is built by the fitted helix builder (not available yet)` (TW-07).

## Wave 3: W3-A and W3-B done; W3-C remaining

- W3-A, the density fit in the window: *Global search* and *Refine* on the Structure page run the
  command line's fit as a window job, the placement stored as one undo step, nothing on Cancel, a map
  changed during the fit not given the old placement, Render / Render still keep an automatic fit's
  placement (DIFFERENCES WF-01 ... WF-03; GUI case `canvas/density_fit`).
- W3-B, the PDB clean-up and structure builder windows behind the host interface (D1), opened from
  the Workflow menu (DIFFERENCES TW-01 ... TW-08; GUI cases `tools/*`; ctest `gui.tools`,
  `gui.tools.scale150`, `gui.tools.scale200`). The non-WC builds of the builder window wait for W2-A
  (`ToolTests::builderNonWcHelixBuildsInTheBackground` is marked `PENDING(W2-A)`).
- Integrated together: a GUI case's `fitted` (`{file: tolerance}`, need `fitted`) is one mechanism
  for both (tests/data/README.md "gui"): fitted files line by line with their decimals within the last
  printed digit or the tolerance, a project's `[density]` table within D10's angle / shift / CC
  tolerances, and the case's log, snapshots and reports with their decimals within the last digit.
- W3-C, the docs (DIFFERENCES "Chemistry" with every CH item, README / BUILDING, LICENSING and
  ATTRIBUTION rows, ARCHITECTURE status), remains; it closes the milestone after W2-A.
