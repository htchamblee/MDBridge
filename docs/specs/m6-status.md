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

## Wave 1: in progress

W1-A PDB hygiene, W1-C dictionaries, W1-D Watson-Crick helix, W1-E density fit; W1-B annotation
after W1-A.
