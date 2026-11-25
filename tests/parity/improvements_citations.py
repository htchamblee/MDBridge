"""MDBridge's citation improvements over MDATK v0.0.4.1 (docs/DIFFERENCES.md, "Citations"), as
changes of MDATK's citation tables and a replacement of its pipeline_resources, so that the
generated table (tools/gen_citation_data.py -> src/core/citation_data.inc) and the expected data of
the generators that record citations (make_citation_cases.py, cli_parity.py, make_render_cli_cases.py
through render_common.setup) are MDBridge's. Standard library only; MDATK is imported by
apply_mdatk() alone.

  CIT-01  MDATK cites MDAnalysis, NumPy, SciPy and Matplotlib for every project because it ran on
          them. MDBridge does not use them, so they are not cited (owner decision: cite the
          originators of the methods MDBridge uses, not one implementation of them). Their
          resources and references leave the table. Methods with a published originator that
          determine a reported result are cited instead, only when the project uses them: the QCP
          superposition (Theobald 2005; Liu, Agrafiotis and Theobald 2010) for an RMSD target
          superposed without a fit selection, Kabsch's superposition (Kabsch 1976, 1978) for a
          target with a fit selection.
  CIT-02  Barnaba is not cited: MDBridge never runs it (with analyses.basepairs.engine "barnaba"
          its builtin detector runs; docs/DIFFERENCES.md "No Barnaba"). Its table entry stays.
  CIT-04  The references of the energy analysis' parameters without parm7 (Amber force-field residue
          templates, Amber atom types) join the table where MDATK's lacks them (ADDED_REFERENCES); the
          analysis cites them in its log line and lj_references, by what it used.
  RN-16   (kept here too) VTK is never cited: "vtk" means MDBridge's Vulkan renderer.

Every function mirrors the C++ rule (src/core/citations.cpp pipelineResources) exactly.
"""

from __future__ import annotations

# CIT-01: resources (and their references) MDBridge does not use
REMOVED_RESOURCES = ("mdanalysis", "numpy", "scipy", "matplotlib")
REMOVED_REFERENCES = ("michaudagrawal2011", "gowers2016", "harris2020", "virtanen2020", "hunter2007")

# CIT-01: the originators of the superposition methods of MDBridge's RMSD (appended to the table).
# `verified`: crossref:<url> is the publisher's metadata (title, authors, journal, volume, issue,
# pages, dates) as deposited with Crossref, read on 2026-11-25; doi.org resolves each DOI to the
# publisher's article page (IUCr Journals, Wiley Online Library).
ADDED_REFERENCES = [
    {
        "key": "theobald2005",
        "text": "Theobald D. L. Acta Crystallogr. A 2005, 61, 478-480. doi:10.1107/S0108767305015266",
        "bibtex": (
            "@article{theobald2005,\n  author = {Theobald, D. L.},\n  title = {{Rapid calculation of RMSDs "
            "using a quaternion-based characteristic polynomial}},\n  journal = {Acta Crystallographica "
            "Section A},\n  year = {2005},\n  volume = {61},\n  number = {4},\n  pages = {478--480},\n  "
            "doi = {10.1107/S0108767305015266}\n}"
        ),
        "doi": "10.1107/S0108767305015266",
        "url": "",
        "verified": "crossref:https://api.crossref.org/works/10.1107/S0108767305015266",
    },
    {
        "key": "liu2010",
        "text": (
            "Liu P., Agrafiotis D. K., Theobald D. L. J. Comput. Chem. 2010, 31, 1561-1563. "
            "doi:10.1002/jcc.21439"
        ),
        "bibtex": (
            "@article{liu2010,\n  author = {Liu, P. and Agrafiotis, D. K. and Theobald, D. L.},\n  title = "
            "{{Fast determination of the optimal rotational matrix for macromolecular superpositions}},\n  "
            "journal = {Journal of Computational Chemistry},\n  year = {2010},\n  volume = {31},\n  number = "
            "{7},\n  pages = {1561--1563},\n  doi = {10.1002/jcc.21439}\n}"
        ),
        "doi": "10.1002/jcc.21439",
        "url": "",
        "verified": "crossref:https://api.crossref.org/works/10.1002/jcc.21439",
    },
    {
        "key": "kabsch1976",
        "text": "Kabsch W. Acta Crystallogr. A 1976, 32, 922-923. doi:10.1107/S0567739476001873",
        "bibtex": (
            "@article{kabsch1976,\n  author = {Kabsch, W.},\n  title = {{A solution for the best rotation "
            "to relate two sets of vectors}},\n  journal = {Acta Crystallographica Section A},\n  year = "
            "{1976},\n  volume = {32},\n  number = {5},\n  pages = {922--923},\n  doi = "
            "{10.1107/S0567739476001873}\n}"
        ),
        "doi": "10.1107/S0567739476001873",
        "url": "",
        "verified": "crossref:https://api.crossref.org/works/10.1107/S0567739476001873",
    },
    {
        "key": "kabsch1978",
        "text": "Kabsch W. Acta Crystallogr. A 1978, 34, 827-828. doi:10.1107/S0567739478001680",
        "bibtex": (
            "@article{kabsch1978,\n  author = {Kabsch, W.},\n  title = {{A discussion of the solution for "
            "the best rotation to relate two sets of vectors}},\n  journal = {Acta Crystallographica "
            "Section A},\n  year = {1978},\n  volume = {34},\n  number = {5},\n  pages = {827--828},\n  "
            "doi = {10.1107/S0567739478001680}\n}"
        ),
        "doi": "10.1107/S0567739478001680",
        "url": "",
        "verified": "crossref:https://api.crossref.org/works/10.1107/S0567739478001680",
    },
    # CIT-04: the force fields of the energy analysis' parameters without parm7 (data/sources/forcefields,
    # tools/gen_energy_lj.py): the papers the Amber 2024 Reference Manual gives for ff14SB, OL15 and
    # OL3 and for the Li-Merz IOD ion set, as far as MDATK's table does not have them already
    # (cornell1995, wang2000, perez2007, zgarbova2011, jorgensen1983, joung2008, li2013). The energy
    # analysis names those it used in its log and in its metadata (lj_references); they are not
    # project resources, since which ones a run uses depends on the residues it meets.
    {
        "key": "maier2015",
        "text": (
            "Maier J. A., Martinez C., Kasavajhala K., Wickstrom L., Hauser K. E., Simmerling C. J. Chem. "
            "Theory Comput. 2015, 11, 3696-3713. doi:10.1021/acs.jctc.5b00255"
        ),
        "bibtex": (
            "@article{maier2015,\n  author = {Maier, J. A. and Martinez, C. and Kasavajhala, K. and "
            "Wickstrom, L. and Hauser, K. E. and Simmerling, C.},\n  title = {{ff14SB: Improving the "
            "Accuracy of Protein Side Chain and Backbone Parameters from ff99SB}},\n  journal = {Journal of "
            "Chemical Theory and Computation},\n  year = {2015},\n  volume = {11},\n  number = {8},\n  "
            "pages = {3696--3713},\n  doi = {10.1021/acs.jctc.5b00255}\n}"
        ),
        "doi": "10.1021/acs.jctc.5b00255",
        "url": "",
        "verified": "crossref:https://api.crossref.org/works/10.1021/acs.jctc.5b00255",
    },
    {
        "key": "banas2010",
        "text": (
            "Banáš P., Hollas D., Zgarbová M., Jurečka P., Orozco M., Cheatham T. E., Šponer J., Otyepka M. "
            "J. Chem. Theory Comput. 2010, 6, 3836-3849. doi:10.1021/ct100481h"
        ),
        "bibtex": (
            "@article{banas2010,\n  author = {Banáš, P. and Hollas, D. and Zgarbová, M. and Jurečka, P. and "
            "Orozco, M. and Cheatham, III, T. E. and Šponer, J. and Otyepka, M.},\n  title = {{Performance "
            "of Molecular Mechanics Force Fields for {RNA} Simulations: Stability of {UUCG} and {GNRA} "
            "Hairpins}},\n  journal = {Journal of Chemical Theory and Computation},\n  year = {2010},\n  "
            "volume = {6},\n  number = {12},\n  pages = {3836--3849},\n  doi = {10.1021/ct100481h}\n}"
        ),
        "doi": "10.1021/ct100481h",
        "url": "",
        "verified": "crossref:https://api.crossref.org/works/10.1021/ct100481h",
    },
    {
        "key": "krepl2012",
        "text": (
            "Krepl M., Zgarbová M., Stadlbauer P., Otyepka M., Banáš P., Koča J., Cheatham T. E., Jurečka "
            "P., Šponer J. J. Chem. Theory Comput. 2012, 8, 2506-2520. doi:10.1021/ct300275s"
        ),
        "bibtex": (
            "@article{krepl2012,\n  author = {Krepl, M. and Zgarbová, M. and Stadlbauer, P. and Otyepka, M. "
            "and Banáš, P. and Koča, J. and Cheatham, III, T. E. and Jurečka, P. and Šponer, J.},\n  title "
            "= {{Reference Simulations of Noncanonical Nucleic Acids with Different $\\chi$ Variants of the "
            "{AMBER} Force Field: Quadruplex {DNA}, Quadruplex {RNA}, and {Z-DNA}}},\n  journal = {Journal "
            "of Chemical Theory and Computation},\n  year = {2012},\n  volume = {8},\n  number = {7},\n  "
            "pages = {2506--2520},\n  doi = {10.1021/ct300275s}\n}"
        ),
        "doi": "10.1021/ct300275s",
        "url": "",
        "verified": "crossref:https://api.crossref.org/works/10.1021/ct300275s",
    },
    {
        "key": "zgarbova2013",
        "text": (
            "Zgarbová M., Luque F. J., Šponer J., Cheatham T. E., Otyepka M., Jurečka P. J. Chem. Theory "
            "Comput. 2013, 9, 2339-2354. doi:10.1021/ct400154j"
        ),
        "bibtex": (
            "@article{zgarbova2013,\n  author = {Zgarbová, M. and Luque, F. J. and Šponer, J. and Cheatham, "
            "III, T. E. and Otyepka, M. and Jurečka, P.},\n  title = {{Toward Improved Description of {DNA} "
            "Backbone: Revisiting Epsilon and Zeta Torsion Force Field Parameters}},\n  journal = {Journal "
            "of Chemical Theory and Computation},\n  year = {2013},\n  volume = {9},\n  number = {5},\n  "
            "pages = {2339--2354},\n  doi = {10.1021/ct400154j}\n}"
        ),
        "doi": "10.1021/ct400154j",
        "url": "",
        "verified": "crossref:https://api.crossref.org/works/10.1021/ct400154j",
    },
    {
        "key": "zgarbova2015",
        "text": (
            "Zgarbová M., Šponer J., Otyepka M., Cheatham T. E., Galindo-Murillo R., Jurečka P. J. Chem. "
            "Theory Comput. 2015, 11, 5723-5736. doi:10.1021/acs.jctc.5b00716"
        ),
        "bibtex": (
            "@article{zgarbova2015,\n  author = {Zgarbová, M. and Šponer, J. and Otyepka, M. and Cheatham, "
            "III, T. E. and Galindo-Murillo, R. and Jurečka, P.},\n  title = {{Refinement of the "
            "Sugar--Phosphate Backbone Torsion Beta for {AMBER} Force Fields Improves the Description of "
            "{Z-} and {B-DNA}}},\n  journal = {Journal of Chemical Theory and Computation},\n  year = "
            "{2015},\n  volume = {11},\n  number = {12},\n  pages = {5723--5736},\n  doi = "
            "{10.1021/acs.jctc.5b00716}\n}"
        ),
        "doi": "10.1021/acs.jctc.5b00716",
        "url": "",
        "verified": "crossref:https://api.crossref.org/works/10.1021/acs.jctc.5b00716",
    },
    {
        "key": "galindomurillo2016",
        "text": (
            "Galindo-Murillo R., Robertson J. C., Zgarbová M., Šponer J., Otyepka M., Jurečka P., Cheatham "
            "T. E. J. Chem. Theory Comput. 2016, 12, 4114-4127. doi:10.1021/acs.jctc.6b00186"
        ),
        "bibtex": (
            "@article{galindomurillo2016,\n  author = {Galindo-Murillo, R. and Robertson, J. C. and "
            "Zgarbová, M. and Šponer, J. and Otyepka, M. and Jurečka, P. and Cheatham, III, T. E.},\n  "
            "title = {{Assessing the Current State of {Amber} Force Field Modifications for {DNA}}},\n  "
            "journal = {Journal of Chemical Theory and Computation},\n  year = {2016},\n  volume = {12},\n  "
            "number = {8},\n  pages = {4114--4127},\n  doi = {10.1021/acs.jctc.6b00186}\n}"
        ),
        "doi": "10.1021/acs.jctc.6b00186",
        "url": "",
        "verified": "crossref:https://api.crossref.org/works/10.1021/acs.jctc.6b00186",
    },
    {
        "key": "li2015",
        "text": "Li P., Song L. F., Merz K. M. J. Phys. Chem. B 2015, 119, 883-895. doi:10.1021/jp505875v",
        "bibtex": (
            "@article{li2015,\n  author = {Li, P. and Song, L. F. and Merz, Jr., K. M.},\n  title = "
            "{{Parameterization of Highly Charged Metal Ions Using the 12-6-4 {LJ}-Type Nonbonded Model in "
            "Explicit Water}},\n  journal = {The Journal of Physical Chemistry B},\n  year = {2015},\n  "
            "volume = {119},\n  number = {3},\n  pages = {883--895},\n  doi = {10.1021/jp505875v}\n}"
        ),
        "doi": "10.1021/jp505875v",
        "url": "",
        "verified": "crossref:https://api.crossref.org/works/10.1021/jp505875v",
    },
]

ADDED_RESOURCES = [
    {
        "id": "qcp",
        "name": "QCP superposition (RMSD)",
        "keys": ["theobald2005", "liu2010"],
        "note": "RMSD targets superposed without a fit selection (analyses.rmsd).",
    },
    {
        "id": "kabsch",
        "name": "Kabsch superposition (RMSD fit selection)",
        "keys": ["kabsch1976", "kabsch1978"],
        "note": "RMSD targets superposed on a fit selection (analyses.rmsd).",
    },
]


def patch_data(references, resources):
    """MDATK's REFERENCE_DATA / RESOURCE_DATA as MDBridge has them: new lists, MDATK's order with
    CIT-01's removals, the additions at the end."""
    refs = [dict(r) for r in references if r["key"] not in REMOVED_REFERENCES]
    res = [dict(r) for r in resources if r["id"] not in REMOVED_RESOURCES]
    refs += [dict(r) for r in ADDED_REFERENCES]
    res += [dict(r) for r in ADDED_RESOURCES]
    return refs, res


def pipeline_resources(project):
    """core::citations::pipelineResources: cpptraj (cleaning.enabled); qcp / kabsch (analyses.rmsd
    enabled with targets and superposition: qcp for a target whose fit_selection is blank, kabsch
    for one that has a fit selection); ffmpeg (output.video.enabled). Never mdanalysis, numpy,
    scipy, matplotlib (CIT-01), barnaba (CIT-02) or vtk (RN-16)."""
    out = []
    if project.cleaning.enabled:
        out.append("cpptraj")
    rmsd = project.analyses.rmsd
    if rmsd.enabled and rmsd.targets and rmsd.superposition:
        if any(not t.fit_selection.strip() for t in rmsd.targets):
            out.append("qcp")
        if any(t.fit_selection.strip() for t in rmsd.targets):
            out.append("kabsch")
    if project.output.video.enabled:
        out.append("ffmpeg")
    return out


def apply_mdatk():
    """Patches MDATK's citation tables (in place, so every module holding them sees the change) and
    pipeline_resources (mdatk.citations and the name mdatk.pipeline imported)."""
    from mdatk import citations as C

    refs, res = patch_data(
        [{"key": r.key, "text": r.text, "bibtex": r.bibtex, "doi": r.doi, "url": r.url, "verified": r.verified}
         for r in C.REFERENCES.values()],
        [{"id": r.id, "name": r.name, "keys": list(r.keys), "note": r.note} for r in C.RESOURCES.values()])
    C.REFERENCES.clear()
    C.REFERENCES.update({r["key"]: C.Reference(**r) for r in refs})
    C.RESOURCES.clear()
    C.RESOURCES.update({r["id"]: C.Resource(r["id"], r["name"], tuple(r["keys"]), r.get("note", ""))
                        for r in res})
    C.pipeline_resources = pipeline_resources  # project_resources calls it through the module
    from mdatk import pipeline

    pipeline.pipeline_resources = pipeline_resources
