#include "gui/tools/pdb_tools_window.hpp"

#include <exception>
#include <string>
#include <vector>

#include <QCheckBox>
#include <QComboBox>
#include <QLineEdit>
#include <QPushButton>
#include <QScrollArea>
#include <QSpinBox>
#include <QSplitter>
#include <QVBoxLayout>

#include "gui/display_path.hpp"
#include "gui/forms/form_context.hpp"
#include "gui/pages/form_page.hpp"
#include "gui/tools/host.hpp"
#include "gui/tools/tool_actions.hpp"

namespace mdbridge::app::tools {

namespace {

using forms::toQString;
using forms::toStdString;

std::vector<forms::FileFilter> pdbFiles() { return {{QStringLiteral("PDB files"), {"*.pdb", "*.ent", "*.pdb1"}}}; }
std::vector<forms::FileFilter> fastaFiles() { return {{QStringLiteral("FASTA files"), {"*.fasta", "*.fa", "*.fas"}}}; }

}  // namespace

PdbToolsWindow::PdbToolsWindow(Host& host, QWidget* parent)
    : ToolWindow(host, QStringLiteral("pdb_tools"), QStringLiteral("PDB clean-up"), parent) {
    QVBoxLayout* column = nullptr;
    QScrollArea* scroll = makeFormScroll(QStringLiteral("form"), QStringLiteral("PDB clean-up options"), &column, splitter());
    QWidget* body = scroll->widget();
    const QString page = toolId();

    // ---- Files
    auto* files = new pages::Section(page, QStringLiteral("files"), QStringLiteral("Files"), body);
    m_input = makePath(QStringLiteral("input"), QStringLiteral("Input PDB"), QStringLiteral("Select a PDB file"),
                       forms::PathMode::OpenFile, pdbFiles());
    files->form()->add(QStringLiteral("Input PDB"), m_input,
                       QStringLiteral("PDB to repair (ChimeraX-damaged lines are fixed automatically)"));
    m_output = makePath(QStringLiteral("output"), QStringLiteral("Output PDB"), QStringLiteral("Output PDB"),
                        forms::PathMode::SaveFile, pdbFiles());
    files->form()->add(QStringLiteral("Output PDB"), m_output);
    column->addWidget(files);

    // ---- Names and records
    auto* names = new pages::Section(page, QStringLiteral("names"), QStringLiteral("Names and records"), body);
    m_scheme = makeCombo(QStringLiteral("scheme"), QStringLiteral("Residue names"),
                         {{"", "keep the names"},
                          {"amber", "Amber (L: LA/LC/LG/LU + 5'/3', Dudek & Trylska)"},
                          {"ccd", "CCD for Coot / Phenix (L: 0A/0C/0G/0U + 0A5/0A3)"},
                          {"rosetta", "Rosetta (L: 0A/0C/0G/0U, all ATOM)"},
                          {"pdb", "PDB deposition (L: 0A..., non-standard as HETATM)"}},
                         "");
    names->form()->add(QStringLiteral("Residue names"), m_scheme,
                       QStringLiteral("Base from the name, chirality from the sugar geometry (or the ranges below), "
                                      "5'/3' variants from the segment ends"));
    m_terminals = makeCheck(QStringLiteral("terminals"), QStringLiteral("5'/3' terminal names (LA5/LA3, A5/A3 …)"), true,
                            QStringLiteral("Off: chain ends keep the internal names (LA, A), for leap libraries that "
                                           "map terminals with addPdbResMap"));
    names->form()->addWide(m_terminals);
    m_detect = makeCheck(QStringLiteral("detect"),
                         QStringLiteral("Detect L-nucleotides from the sugar chirality (C1', C3', C4')"), true);
    names->form()->addWide(m_detect);
    m_lRanges = makeLine(QStringLiteral("l_ranges"), QStringLiteral("Force L for"),
                         QStringLiteral("e.g. C:45-73 or C (whole chain)"));
    names->form()->add(QStringLiteral("Force L for"), m_lRanges,
                       QStringLiteral("Residue ranges that are L whatever the geometry says (overrides the detection)"));
    m_records = makeCombo(QStringLiteral("records"), QStringLiteral("ATOM / HETATM"),
                          {{"", "follow the name scheme"},
                           {"keep", "keep ATOM / HETATM"},
                           {"rosetta", "all polymer residues as ATOM (Rosetta)"},
                           {"pdb", "non-standard residues as HETATM (PDB)"}},
                          "");
    names->form()->add(QStringLiteral("ATOM / HETATM"), m_records);
    m_keepConect = makeCheck(QStringLiteral("keep_conect"),
                             QStringLiteral("Keep CONECT records (stripped until final deposition processing)"));
    names->form()->addWide(m_keepConect);
    m_strip5p = makeCheck(QStringLiteral("strip_5p"), QStringLiteral("Remove 5'-terminal phosphates (for Amber)"), false,
                          QStringLiteral("Drop P/OP1/OP2/OP3 of every chain's first nucleotide: Amber's 5' units (G5, "
                                         "LG5 ...) have a 5'-OH. Hydrogenated residues get HO5'."));
    names->form()->addWide(m_strip5p);
    connect(m_scheme, &QComboBox::currentIndexChanged, this, [this] {
        m_strip5p->setChecked(m_strip5p->isChecked() || currentId(m_scheme) == QLatin1String("amber"));
    });
    column->addWidget(names);

    // ---- Chains
    auto* chains = new pages::Section(page, QStringLiteral("chains"), QStringLiteral("Chains"), body);
    m_chainMode = makeCombo(QStringLiteral("chain_mode"), QStringLiteral("Chain ids"),
                            {{"keep", "keep the chain ids"},
                             {"reference", "copy from a reference PDB"},
                             {"list", "ids per TER segment (list)"},
                             {"assign", "A, B, C … per TER segment"}},
                            "keep");
    chains->form()->add(QStringLiteral("Chain ids"), m_chainMode,
                        QStringLiteral("cpptraj writes PDBs without chain ids: copy them from the structure the "
                                       "topology was built from, give one id per TER segment, or assign A, B, C …"));
    m_reference = makePath(QStringLiteral("reference"), QStringLiteral("Reference PDB"), QStringLiteral("Reference PDB"),
                           forms::PathMode::OpenFile, pdbFiles());
    chains->form()->add(QStringLiteral("Reference PDB"), m_reference);
    m_copyNumbers = makeCheck(QStringLiteral("copy_numbers"), QStringLiteral("Also copy the residue numbers"));
    chains->form()->addWide(m_copyNumbers);
    m_chainIds = makeLine(QStringLiteral("chain_ids"), QStringLiteral("Ids per segment"), QStringLiteral("e.g. A,B,C"));
    chains->form()->add(QStringLiteral("Ids per segment"), m_chainIds);
    m_splitSegments = makeCheck(QStringLiteral("split_segments"), QStringLiteral("Give every TER segment its own new chain id"));
    chains->form()->addWide(m_splitSegments);
    m_splitAt = makeLine(QStringLiteral("split_at"), QStringLiteral("Split before"), QStringLiteral("e.g. A:31"));
    chains->form()->add(QStringLiteral("Split before"), m_splitAt,
                        QStringLiteral("Start a new chain (new id) at these residues"));
    m_combine = makeLine(QStringLiteral("combine"), QStringLiteral("Combine chains"), QStringLiteral("e.g. A,B"));
    m_combineInto = makeLine(QStringLiteral("combine_into"), QStringLiteral("Combine into"), QStringLiteral("first"),
                             QString(), 4);
    chains->form()->addRow(QStringLiteral("Combine chains / into"), {m_combine, m_combineInto});
    m_sortChains = makeCheck(QStringLiteral("sort_chains"), QStringLiteral("Write the chains in alphabetical order"));
    chains->form()->addWide(m_sortChains);
    connect(m_chainMode, &QComboBox::currentIndexChanged, this, &PdbToolsWindow::sync);
    column->addWidget(chains);

    // ---- Numbering and sequence
    auto* numbering = new pages::Section(page, QStringLiteral("numbering"), QStringLiteral("Numbering and sequence"), body);
    m_renumber = makeCheck(QStringLiteral("renumber"), QStringLiteral("Renumber residues (continuous across chains)"));
    numbering->form()->addWide(m_renumber);
    m_start = makeSpin(QStringLiteral("start"), QStringLiteral("First number"), -999, 9999, 1);
    numbering->form()->add(QStringLiteral("First number"), m_start);
    m_perChain = makeCheck(QStringLiteral("per_chain"), QStringLiteral("Restart the numbering in every chain"));
    numbering->form()->addWide(m_perChain);
    m_keepGaps = makeCheck(QStringLiteral("keep_gaps"), QStringLiteral("Keep numbering gaps"));
    numbering->form()->addWide(m_keepGaps);
    m_seqres = makeCheck(QStringLiteral("seqres"), QStringLiteral("Write SEQRES records"));
    numbering->form()->addWide(m_seqres);
    m_fasta = makePath(QStringLiteral("fasta"), QStringLiteral("FASTA (optional)"), QStringLiteral("FASTA file"),
                       forms::PathMode::SaveFile, fastaFiles());
    numbering->form()->add(QStringLiteral("FASTA (optional)"), m_fasta);
    connect(m_renumber, &QCheckBox::toggled, this, &PdbToolsWindow::sync);
    column->addWidget(numbering);

    // ---- Annotation (CH-14: MDBridge has no Barnaba; every engine is the builtin detector)
    auto* annotation = new pages::Section(page, QStringLiteral("annotation"), QStringLiteral("Annotation"), body);
    m_annotate = makeCheck(QStringLiteral("annotate"), QStringLiteral("Annotate base pairs as REMARK 950 records"), false,
                           QStringLiteral("Base pairs, stacking and dot-bracket written into the header of the output "
                                          "(as 'mdbridge pdb annotate')"));
    annotation->form()->addWide(m_annotate);
    m_engine = makeCombo(QStringLiteral("ann_engine"), QStringLiteral("Engine"),
                         {{"auto", "automatic: the builtin detector"},
                          {"builtin", "builtin (canonical pairs only)"},
                          {"barnaba", "Barnaba - not available (the builtin detector, with a warning)"}},
                         "auto");
    annotation->form()->add(QStringLiteral("Engine"), m_engine,
                            QStringLiteral("MDBridge has no Barnaba: every engine is the builtin detector, which finds "
                                           "Watson-Crick and G-U wobble pairs from H-bond distances"));
    m_annPairs = makeCheck(QStringLiteral("ann_pairs"), QStringLiteral("Base pairs (PAIR)"), true);
    m_annStacking = makeCheck(QStringLiteral("ann_stacking"), QStringLiteral("Stacking (STACK, Barnaba only)"), true,
                              QStringLiteral("Barnaba's stacking records: MDBridge writes none (a note says so)"));
    m_annDbn = makeCheck(QStringLiteral("ann_dbn"), QStringLiteral("Sequence and dot-bracket (SEQ / DBN)"), true);
    for (QCheckBox* box : {m_annPairs, m_annStacking, m_annDbn}) annotation->form()->addWide(box);
    connect(m_annotate, &QCheckBox::toggled, this, &PdbToolsWindow::sync);
    column->addWidget(annotation);
    column->addStretch(1);
    splitter()->addWidget(scroll);

    addReport();
    addButton(QStringLiteral("check"), QStringLiteral("Check input"),
              QStringLiteral("Chains, sequences, residue names and the sugar chirality of every nucleotide (nothing is "
                             "written)"));
    QPushButton* write = addButton(QStringLiteral("write"), QStringLiteral("Write cleaned PDB"),
                                   QStringLiteral("Apply the options and write the output file"));
    write->setDefault(true);  // GTK's suggested-action
    QPushButton* close = addButton(QStringLiteral("close"), QStringLiteral("Close"));
    connect(qobject_cast<QPushButton*>(control(QStringLiteral("check"))), &QPushButton::clicked, this,
            &PdbToolsWindow::runCheck);
    connect(write, &QPushButton::clicked, this, &PdbToolsWindow::runClean);
    connect(close, &QPushButton::clicked, this, &QWidget::close);
    connect(m_input, &PathEntry::pathChanged, this, &PdbToolsWindow::onInput);
    finishLayout(980, 720);
    sync();
}

void PdbToolsWindow::sync() {
    const QString mode = currentId(m_chainMode);
    m_reference->setEnabled(mode == QLatin1String("reference"));
    m_copyNumbers->setEnabled(mode == QLatin1String("reference"));
    m_chainIds->setEnabled(mode == QLatin1String("list"));
    for (QWidget* w : std::initializer_list<QWidget*>{m_start, m_perChain, m_keepGaps}) w->setEnabled(m_renumber->isChecked());
    for (QWidget* w : std::initializer_list<QWidget*>{m_engine, m_annPairs, m_annStacking, m_annDbn}) {
        w->setEnabled(m_annotate->isChecked());
    }
}

// <input>_clean.pdb while the output is empty or still the previous suggestion
void PdbToolsWindow::onInput() {
    const QString current = m_output->path();
    if (!current.isEmpty() && current != m_autoOutput) {
        return;  // chosen by hand
    }
    m_autoOutput = toQString(suggestedOutput(toStdString(m_input->path()), "clean"));
    m_output->setPath(m_autoOutput);
}

PdbCleanRequest PdbToolsWindow::request(const QString& input, const QString& output) const {
    PdbCleanRequest r;
    r.input = toStdString(input);
    r.output = toStdString(output);
    r.fasta = toStdString(m_fasta->path());
    const QString mode = currentId(m_chainMode);
    r.scheme = toStdString(currentId(m_scheme));
    r.records = toStdString(currentId(m_records));
    r.stripConect = !m_keepConect->isChecked();
    r.detectL = m_detect->isChecked();
    r.lRanges = toStdString(m_lRanges->text().trimmed());
    if (mode == QLatin1String("reference")) r.chainReference = toStdString(m_reference->path());
    if (mode == QLatin1String("list")) r.chainIds = splitIds(toStdString(m_chainIds->text()));
    r.assignChains = mode == QLatin1String("assign");
    r.copyNumbers = mode == QLatin1String("reference") && m_copyNumbers->isChecked();
    r.splitSegments = m_splitSegments->isChecked();
    r.splitAt = toStdString(m_splitAt->text().trimmed());
    r.combine = splitIds(toStdString(m_combine->text())).value_or(std::vector<std::string>{});
    if (const QString into = m_combineInto->text().trimmed(); !into.isEmpty()) r.combineInto = toStdString(into);
    r.sortChains = m_sortChains->isChecked();
    r.renumber = m_renumber->isChecked();
    r.start = m_start->value();
    r.perChain = m_perChain->isChecked();
    r.keepGaps = m_keepGaps->isChecked();
    r.seqres = m_seqres->isChecked();
    r.terminalNames = m_terminals->isChecked();
    r.strip5p = m_strip5p->isChecked();
    if (m_annotate->isChecked()) {
        chem::AnnotateOptions a;
        a.pairs = m_annPairs->isChecked();
        a.stacking = m_annStacking->isChecked();
        a.dotbracket = m_annDbn->isChecked();
        a.engine = toStdString(currentId(m_engine));
        if (a.engine.empty()) a.engine = "auto";
        r.annotation = a;
    }
    return r;
}

QString PdbToolsWindow::inputFile() {
    const QString path = m_input->path();
    if (path.isEmpty() || !isFile(toStdString(path))) {
        showReport(QStringLiteral("Choose an existing input PDB file first."));
        return {};
    }
    return path;
}

void PdbToolsWindow::runCheck() {
    const QString path = inputFile();
    if (path.isEmpty()) {
        return;
    }
    try {
        showReport(toQString(pdbCheckReport(toStdString(path), toStdString(m_lRanges->text().trimmed()))));
    } catch (const std::exception& e) {  // TW-02 (MDATK 0.1.0): also a file that is not UTF-8 text
        showReport(QStringLiteral("ERROR: ") + toQString(errorText(e)));
    }
}

void PdbToolsWindow::runClean() {
    const QString path = inputFile();
    if (path.isEmpty()) {
        return;
    }
    const QString out = m_output->path();
    if (out.isEmpty()) {
        showReport(QStringLiteral("Choose an output file."));
        return;
    }
    if (samePath(toStdString(out), toStdString(path))) {
        showReport(QStringLiteral("The output must not overwrite the input file."));
        return;
    }
    ToolResult result;
    try {
        result = cleanPdb(request(path, out), host().projectPath());
    } catch (const std::exception& e) {
        showReport(QStringLiteral("ERROR: ") + toQString(errorText(e)));
        return;
    }
    QStringList lines{QStringLiteral("Wrote ") + gui::displayPath(out)};
    if (!m_fasta->path().isEmpty()) lines << QStringLiteral("Wrote ") + gui::displayPath(m_fasta->path());
    lines << QString();
    for (const std::string& note : result.notes) lines << toQString(note);
    showReport(lines.join(QLatin1Char('\n')));
    host().appendLog(QStringLiteral("PDB clean-up: wrote ") + gui::displayPath(out));
}

}  // namespace mdbridge::app::tools
