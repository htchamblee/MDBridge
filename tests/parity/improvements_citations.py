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
  CIT-04  The energy analysis' parameters without parm7 (Amber force-field residue templates, Amber
          atom types) cite the originators of the numbers used; the one reference MDATK's table lacks
          (Li et al. 2015) joins it; the analysis cites them in its log line and lj_references. MDATK
          0.1.1's resource generic_lj is known (the typed layer's references) but never a project
          resource.
  CIT-06  Run-time citations: the references a run's cached results record (meta "citation_keys":
          the energy analysis' force fields and atom types, Cordero et al.'s radii of a bond guess)
          follow the settings' ones in the parameter log and the render step's log entry
          (_patch_run_citations); Cordero et al. 2008 joins the table.
  RN-16   (kept here too) VTK is never cited: "vtk" means MDBridge's Vulkan renderer.

Every function mirrors the C++ rule (src/core/citations.cpp pipelineResources) exactly.
"""

from __future__ import annotations

import os

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
    # CIT-04: the energy analysis' parameters without parm7 (data/sources/forcefields,
    # tools/gen_energy_lj.py) cite who originated the numbers used: Cornell et al. 1995 (the
    # charges and Lennard-Jones parameters of the protein and nucleic templates and of parm94's
    # types), Wang et al. 2000 (Cl, Br), Jorgensen et al. 1983 (TIP3P), Joung & Cheatham 2008 and
    # Li et al. 2013 (ions), all in MDATK's table already, and Li et al. 2015 (the 3+ / 4+ ions of
    # tip3p.xml), added here. The analysis names those it used in its log and in its metadata
    # (lj_references); they are no project resource, since which ones a run uses depends on the
    # residues it meets.
    # CIT-06: the covalent radii of the analyses' bond guess (data/sources/covalent_radii), cited by the
    # runs that guess bonds (the energy and atoms analyses' citation_keys; CIT-03, CIT-05)
    {
        "key": "cordero2008",
        "text": (
            "Cordero B., Gómez V., Platero-Prats A. E., Revés M., Echeverría J., Cremades E., Barragán F., "
            "Alvarez S. Dalton Trans. 2008, 2832-2838. doi:10.1039/b801115j"
        ),
        "bibtex": (
            "@article{cordero2008,\n  author = {Cordero, B. and Gómez, V. and Platero-Prats, A. E. and Revés, M. "
            "and Echeverría, J. and Cremades, E. and Barragán, F. and Alvarez, S.},\n  title = {{Covalent radii "
            "revisited}},\n  journal = {Dalton Transactions},\n  year = {2008},\n  number = {21},\n  pages = "
            "{2832--2838},\n  doi = {10.1039/b801115j}\n}"
        ),
        "doi": "10.1039/b801115j",
        "url": "",
        "verified": "crossref:https://api.crossref.org/works/10.1039/b801115j",
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
    # CIT-04: MDATK 0.1.1's resource for its generic Lennard-Jones table (citations.py, generic_lj.py
    # CITATION_KEYS), known here so that its params logs and project citations that name it resolve.
    # MDBridge cites by run (lj_references) and never lists it in pipeline_resources; its keys are those
    # of MDBridge's typed layer (docs/GENERIC_LJ.md): MDATK's four plus Jorgensen et al. 1983 (TIP3P's
    # water oxygen, type OW).
    {
        "id": "generic_lj",
        "name": "Generic Lennard-Jones table of the energy analysis (Amber atom types and ions)",
        "keys": ["cornell1995", "wang2000", "jorgensen1983", "joung2008", "li2013"],
        "note": "topologies without atom types; docs/GENERIC_LJ.md",
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
    _patch_run_citations(pipeline, C)


def with_run_keys(settings_keys, run_keys, references):
    """CIT-06: core::citations::withRunKeys: the settings' keys, then the run's, without repeats and
    keys the table lacks."""
    out = []
    for key in list(settings_keys) + list(run_keys):
        if key in references and key not in out:
            out.append(key)
    return out


def run_citation_keys(pipe):
    """CIT-06: analysis::Pipeline::runCitationKeys: the "citation_keys" of the cached results of the
    last analyze() (analysis_keys order), without repeats."""
    out = []
    for name, key in pipe.analysis_keys.items():
        loaded = pipe.cache.load_result(name, key)
        if loaded is None:
            continue
        for k in loaded[0].meta.get("citation_keys", []) or []:
            if k not in out:
                out.append(k)
    return out


def _patch_run_citations(pipeline, C):
    """CIT-06: the parameter log and the render step's log entry cite what the run's cached results
    record after the settings' citations (MDATK's write_params_log and run with that one change)."""
    if getattr(pipeline.Pipeline, "_cit06", False):
        return
    from mdatk.projectfile import logged

    def write_params_log(self, outputs):
        from mdatk.projectfile import write_params_log as write

        path = os.path.join(self.project.out_dir, pipeline.PARAMS_LOG_NAME)
        inputs = self._inputs
        run = {
            "project": pipeline.PROJECT_TITLE,
            "environment": pipeline.environment_info(self.project),
            "inputs": {
                "topology": inputs.topology if inputs else None,
                "trajectories": inputs.trajectories if inputs else None,
                "cleaning": ({"key": inputs.cleaned.key, "script": inputs.cleaned.script,
                              "log": inputs.cleaned.log, "reused": inputs.cleaned.reused}
                             if inputs and inputs.cleaned else None),
            },
            "cache_dir": self.project.cache_dir,
            "analysis_keys": self.analysis_keys,
            "timestamp": self._stamp,
            "outputs": outputs.summary(),
        }
        keys = with_run_keys(C.keys_for(C.project_resources(self.project)), run_citation_keys(self), C.REFERENCES)
        run["citation_keys"] = keys
        run["citations"] = C.text_lines(keys)
        write(path, self.project, run)
        return path

    def run(self, video=None, charts=None, pdb=None, force_analysis=False):
        p = self.project
        video = p.output.video.enabled if video is None else video
        charts = p.output.charts.enabled if charts is None else charts
        pdb = p.output.pdb.enabled if pdb is None else pdb
        self.check()
        with logged(self.project_path, "render", "Analyses and outputs",
                    C.keys_for(pipeline.pipeline_resources(p))) as entry:
            outputs = self._run(video, charts, pdb, force_analysis)
            entry.citations = with_run_keys(entry.citations, run_citation_keys(self), C.REFERENCES)
            entry.files = [f for f in [outputs.video, *outputs.charts, *outputs.pdbs, outputs.params_log] if f]
            entry.summary = (f"{outputs.dataset.n_frames} frames analysed" if outputs.dataset else "no analyses") + \
                f"; {len(entry.files)} file(s) written to {p.out_dir}"
        return outputs

    pipeline.Pipeline.write_params_log = write_params_log
    pipeline.Pipeline.run = run
    pipeline.Pipeline._cit06 = True
