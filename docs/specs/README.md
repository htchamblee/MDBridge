# Clean-room behaviour specifications

MDBridge reproduces the behaviour of MDAnalysis (LGPL), but none of its code may be copied or
translated (`LICENSING.md`). Where MDBridge must match MDAnalysis, the work is split:

1. A **spec writer** reads MDAnalysis and runs it, and writes a behaviour-only description here:
   prose, maths, verbatim user-visible messages and input → output tables observed by running
   MDAnalysis. No code or pseudo-code, no MDAnalysis identifiers outside quoted messages.
2. A separate **implementer** writes MDBridge's code from the spec, the rest of MDBridge and the
   test data only. Implementers never open MDAnalysis, its documentation beyond the spec, or the
   code a rewrite replaces.

These files are the record of that process and the reference for later changes to the code they
describe. They are kept as written by the spec writers; where MDBridge deliberately differs from
MDAnalysis (for example its independent data tables, `data/sources/`), the code and
`tests/data/README.md` say so.

| Files | Milestone | Describes |
|---|---|---|
| `m2-selection-*.md` | 2 | the selection language: value tokens, ranges, `prop`, parser, grid-search rule |
| `m2-parm7-topology.md`, `m2-pdb-*.md`, `m2-gro.md` | 2 | topology readers, PDB atom-name placement |
| `m2-dcd-unitcell.md`, `m2-netcdf-amber.md`, `m2-amber-ascii.md`, `m2-triclinic-box.md` | 2 | trajectory readers and box conversion |
| `m3-plan.md` | 3 | the analysis milestone's plan (inventory of MDATK's analyses and work breakdown) |
| `m3-capped-distance.md` | 3 | the periodic pair search (`capped_distance`: method choice, distances, pair order, tolerances) |
| `m3-rmsd.md` | 3 | the optimal-superposition RMSD (`rms.rmsd`: centring, weights, degenerate inputs, errors) |
| `m3-universe-notes.md` | 3 | trajectory iteration by index arrays, cells, charges, chain IDs, bonds and fragments as the analyses use them |
| `m4-plan.md` | 4 | the rendering milestone's plan: inventory of MDATK's rendering, dependencies and licences, parity strategy (contents only, owner decision D10), the RN-xx improvements with their status, work breakdown and owner decisions. Rendering ports MDATK's own code and permissively licensed Matplotlib / NumPy / SciPy algorithms (with notices); no MDAnalysis behaviour is involved, so it has no clean-room specs |
| `m5-plan.md` | 5 | the interface milestone's plan: inventory of MDATK's analysis window, canvas layout editor, dialogs and the commands still missing, the boundary to milestones 6 and 7, parity strategy for the GUI (project files written by GUI actions, Qt Test), the AP-xx improvements (checked against MDATK's later versions), what MDATK added after v0.0.4.1, work breakdown and owner decisions. No MDAnalysis behaviour is involved, so it has no clean-room specs |
