"""MDBridge's analysis versions, given to MDATK's analyses before the generators record anything.

An analysis' version is part of its cache key (docs/ARCHITECTURE.md, "Analysis cache"); MDBridge
bumps it whenever an improvement can change the analysis' results for real inputs. The expected
data (cache keys in pipeline_cases.json, cli_cases.json, the analysis cases) therefore carries
MDBridge's versions, not MDATK's. This table is the one place on the test side that says what they
are; every improvement module calls apply_mdatk() instead of bumping a version itself, so the
versions do not depend on which modules a generator applies, or in which order.

check_cpp() compares the table with the AnalysisDescriptor of each analysis in src/analysis/*.cpp
and raises on any difference (regenerate_all.py and the generators run it): a bump on one side only
fails the regeneration.

Policy: MDBridge has not been released, so an analysis is bumped once, MDATK's version + 1, for
all the improvements of this release that can change its results; a later improvement in the same
release adds its item to the list below instead of bumping again.

Per analysis, the improvements that bumped it (docs/DIFFERENCES.md; "periodic": AN-01 .. AN-04):
  energy 2     AN-07, AN-08 (exclusions); RD-P12 (float64 parm7 charges); periodic, AN-05 (net
               force); SEL-13 (generic LJ elements); AN-25 - AN-27 (parameters and charges without
               parm7: Amber residue templates, Amber atom types)
  hbonds 2     AN-09, AN-10 (guessed hydrogens); periodic; SEL-13 (donor / acceptor elements);
               distance-only pairs counted once (MDATK 0.1.0)
  rmsd 2       AN-16 .. AN-18; AN-20 (NaN for a blown-up fit frame, MDATK 0.1.2)
  proximity 2  AN-15; periodic; SEL-13 (heavy atoms)
  atoms 3      AN-13, AN-14 (with MDATK 0.1.2's residue rule); AN-20 (blown-up frames unaligned);
               SEL-13 (elements, heavy atoms, bond guess) (MDATK: 2)
  contacts 2   periodic; SEL-13 (heavy_only); a group against itself (MDATK 0.1.0), whose strand
               ends follow TOP-01's chains
  mg 2         periodic; SEL-13 (ligand / occupancy elements)
  basepairs 2  periodic, AN-11; SEL-13 (heavy atoms); NA-01 (L-DNA names, CCD atom names, backbone)
  stacking 2   periodic, AN-12 (ring atoms by name: SEL-13 does not change it); NA-01 (CCD atom names)
  trace 5      SEL-13 (heavy-atom centroids of residues without a trace atom); TOP-01 (chains from
               bonds only when they hold the residues together) (MDATK: 4); its key covers every
               group (MDATK 0.1.0; no new results); NA-01 (the per-row array backbone); AN-20
               (blown-up frames unaligned, MDATK 0.1.2)
"""

import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))

# analysis name -> (MDATK module, class, MDBridge version)
VERSIONS = {
    "basepairs": ("mdatk.analysis.basepairs", "BasePairAnalysis", 2),
    "contacts": ("mdatk.analysis.contacts", "ContactAnalysis", 2),
    "energy": ("mdatk.analysis.energy", "EnergyAnalysis", 2),
    "hbonds": ("mdatk.analysis.hbonds", "HBondAnalysis", 2),
    "mg": ("mdatk.analysis.mg", "MgAnalysis", 2),
    "rmsd": ("mdatk.analysis.rmsd", "RMSDAnalysis", 2),
    "stacking": ("mdatk.analysis.stacking", "StackingAnalysis", 2),
    "trace": ("mdatk.analysis.structure", "TraceAnalysis", 5),
    "proximity": ("mdatk.analysis.structure", "ProximityAnalysis", 2),
    "atoms": ("mdatk.analysis.structure", "AtomsAnalysis", 3),
}

_DESCRIPTOR = re.compile(r'AnalysisDescriptor\(\s*\{\s*"(\w+)"\s*,\s*"[^"]*"\s*,\s*(\d+)\s*,')


def cpp_versions(src=None):
    """{name: version} from the AnalysisDescriptor calls in src/analysis/*.cpp."""
    src = src or os.path.join(ROOT, "src", "analysis")
    found = {}
    for name in sorted(os.listdir(src)):
        if name.endswith(".cpp"):
            with open(os.path.join(src, name), encoding="utf-8") as fh:
                for m in _DESCRIPTOR.finditer(fh.read()):
                    found[m.group(1)] = int(m.group(2))
    return found


def check_cpp(root=None):
    """Raise when the table and the C++ descriptors under root (default: this checkout) disagree.
    Skipped where the sources are absent: regenerate_all.py's fresh tree has no C++ sources, so
    regenerate_all.py runs this check on the checkout before the generators."""
    src = os.path.join(root or ROOT, "src", "analysis")
    cpp = cpp_versions(src) if os.path.isdir(src) else {}
    if not cpp:
        return
    table = {name: v[2] for name, v in VERSIONS.items()}
    if cpp != table:
        raise RuntimeError(f"tests/parity/analysis_versions.py {table} does not match the analysis "
                           f"versions in src/analysis/*.cpp {cpp}")


def apply_mdatk():
    """Give MDATK's analysis classes MDBridge's versions (idempotent; a no-op without MDATK)."""
    import importlib

    check_cpp()
    for module_name, class_name, version in VERSIONS.values():
        try:
            module = importlib.import_module(module_name)
        except ImportError:  # a generator without MDATK
            return
        getattr(module, class_name).version = version
