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

## Wave 1: done (CI green at 8d75a4a)

PDB hygiene (`chem/pdbtools.hpp`), annotation (`chem/annotate.hpp`), dictionaries and SMILES
(`chem/dictionaries.hpp`, `chem/smiles.hpp`), the Watson-Crick helix builder (`chem/helix.hpp`) and
the density fit (`chem/density_fit.hpp`); each matches its recorded MDATK cases (byte for byte or
within D10). The regenerate job has a kernel probe: on runners whose pinned kernels give other bits,
MDATK's path-dependent chem results (density-fit evaluation counts, non-WC closures) are compared
structurally (interim; tighten once W2-A makes the non-WC builds well-conditioned).

## Wave 2: in progress

W2-B mutations and fit-to-map, W2-C commands and pipeline. W2-A (other pairings and triplexes) waits
for the conditioning investigation and the owner's decision.
