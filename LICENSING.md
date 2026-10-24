# Licensing plan (not yet in effect)

MDBridge is all rights reserved until its licence is chosen (see `LICENSE`). MDATK, the Python
predecessor, stays GPL v3; its own code (same copyright holder) may be ported into MDBridge
freely, but no third-party GPL code may be copied or bundled.

## Intended licensing

Not adopted yet; it takes effect only after legal review and a formal change of `LICENSE`.

- **Non-commercial use:** free, under a standard source-available non-commercial licence
  (candidate: PolyForm Noncommercial 1.0.0) with a notice prohibiting biological-weapons use.
- **Commercial use:** a paid commercial licence (draft held by the author; to be reviewed by a
  lawyer). Licensor to be the author's LLC once formed, with the copyright assigned to it.

## Rules until then

1. **Copyright:** all code must be owned by the author. Outside contributions require a
   contributor licence agreement that permits relicensing.
2. **Resources:** no third-party or institutional resources (equipment, networks, funded time)
   go into MDBridge.
3. **Dependencies bundled in binaries must allow proprietary use** (permissive or LGPL):

   | Dependency | Licence | Rule | Used since |
   |---|---|---|---|
   | Qt 6 | LGPL v3 | link dynamically; ship the licence and a source pointer | milestone 1 |
   | toml++ | MIT | OK | milestone 1 |
   | CLI11 | BSD-3-Clause | OK | milestone 1 |
   | Vulkan-Headers | Apache 2.0 / MIT | OK (headers only; the Vulkan loader is the user's) | milestone 1 |
   | chemfiles | BSD-3 | OK | planned |
   | Eigen | MPL 2.0 | OK | planned |
   | nanoflann | BSD-2 | OK | planned |
   | FFmpeg | LGPL (only) | no GPL components (no x264); use OpenH264 / AV1 | planned |
   | Lua | MIT | OK | planned |
   | Python, NumPy, SciPy (interim, embedded) | PSF / BSD-style | OK | planned |
   | MDAnalysis (interim) | check the bundled version's licence | verify before bundling | planned |
   | Barnaba | GPL v3 | never bundled; replaced by MDBridge's own base-pair annotation | - |

   Test-only tools (GoogleTest, BSD-3) are not bundled. External programs MDBridge writes scripts
   for (cpptraj, Amber, GROMACS, Rosetta, Phenix, Coot, ...) are not distributed with it and keep
   their own licences.
