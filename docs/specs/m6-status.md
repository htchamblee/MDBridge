# Milestone 6: status

Plan: `docs/specs/m6-plan.md`; owner decisions D1-D14 accepted (2026-11-26), D15 (CHD-01: guanine
amino hydrogens follow the wwPDB CCD) and D16 (CHD-02: non-WC builds in the fitted 4JRD frame).

## Wave 0: done (CI green at a7cacc5)

Chem generators and `improvements_chem.py` (128 command-line cases pending until the commands
exist), data sources (CCD 0A/0C/0G/0U, PDB 4JRD, Arnott A-RNA; `gen_chem_tables.py`), numerics
library (`src/numerics`: ndimage and densityCentroid moved from render, vendored ducc0 FFT,
Powell/Brent, least_squares TRF with LSMR, pair search; NU-01/02). Generators pin numpy/OpenBLAS
kernels; numerics tests are libm-aware (bit-exact with glibc).

Known fragile result for W2-A: MDATK's non-WC / triplex least-squares closures are ill-conditioned
(~0.8 Å spread from the BLAS kernel choice alone, far beyond D10's tolerance). Treat it as a
fragile result to fix, not to copy; options go to the owner before W2-A finishes.

## Wave 1: done (pending CI)

Integrated on `m6-w1-integ` (merges of W1-A/W1-B, W1-C, W1-D, W1-E): PDB hygiene
(`chem::pdbtools`: CH-01, CH-03..07, CH-13, CH-15), base-pair annotation (`chem::annotate`: CH-14,
D13), restraint dictionaries and SMILES (`chem::dictionaries`, `chem::smiles`: DIC-01, DIC-02), the
Watson-Crick helix builder (`chem::buildHelix`: CH-01, CH-02, CHD-01, D13) and the rigid-body
density fit (`chem::FitProblem`, `fitModel`: CH-10, CH-11, D2; several threads). At the merge the
builder's temporary adapter (src/chem/pdb_names) gave way to W1-A's `chem::convertNames` /
`chem::setRecords`, and its own CH-01 check to `io::writePdb`'s (the same message); `chem::kSchemes`
is defined once (chem/pdbtools.hpp). The command lines (`pdb`, `build`, `dict`, `density fit /
simulate`) stay pending for W2-C.

## Wave 2: next

W2-A other pairings and triplexes (the fragile closures above first), W2-B mutations and fit-to-map,
W2-C commands and pipeline.
