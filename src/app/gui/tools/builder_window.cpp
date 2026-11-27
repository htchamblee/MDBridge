#include "gui/tools/builder_window.hpp"

#include <cstdio>
#include <exception>
#include <filesystem>
#include <initializer_list>
#include <string>
#include <utility>
#include <vector>

#include <QCheckBox>
#include <QComboBox>
#include <QDoubleSpinBox>
#include <QLabel>
#include <QLineEdit>
#include <QProgressBar>
#include <QPushButton>
#include <QScrollArea>
#include <QSpinBox>
#include <QSplitter>
#include <QTabBar>
#include <QTabWidget>
#include <QVBoxLayout>

#include "gui/display_path.hpp"
#include "gui/forms/form_context.hpp"
#include "gui/jobs/job_runner.hpp"
#include "gui/pages/form_page.hpp"
#include "gui/tools/host.hpp"
#include "gui/tools/tool_actions.hpp"
#include "mdbridge/core/paths.hpp"

namespace mdbridge::app::tools {

namespace fs = std::filesystem;

namespace {

using forms::toQString;
using forms::toStdString;

constexpr int kProgressSteps = 1000;
const char* const kPages[] = {"helix", "mutate", "fit"};  // MDATK's PAGES

std::vector<forms::FileFilter> pdbFiles() { return {{QStringLiteral("PDB files"), {"*.pdb", "*.ent"}}}; }
std::vector<forms::FileFilter> mapFiles() { return {{QStringLiteral("Density maps"), {"*.mrc", "*.map", "*.ccp4"}}}; }

// what MDATK's windows show of a failed job: str(error)
QString failureText(const jobs::JobFailure& failure) {
    if (failure.error) {
        try {
            std::rethrow_exception(failure.error);
        } catch (const std::exception& e) {
            return toQString(errorText(e));
        } catch (...) {
        }
    }
    return failure.message;
}

QString basename(const QString& path) { return toQString(core::pathBasename(toStdString(path))); }

}  // namespace

BuilderWindow::BuilderWindow(Host& host, QWidget* parent)
    : ToolWindow(host, QStringLiteral("builder"), QStringLiteral("Structure builder"), parent) {
    m_jobs = new jobs::JobRunner({}, this);  // MDATK's worker: no log lines (Worker(job, lambda _t: None, ...))
    m_tabs = new QTabWidget(splitter());
    m_tabs->setObjectName(nameOf(QStringLiteral("pages")));
    m_tabs->setMinimumWidth(500);  // MDATK: notebook.set_size_request(500, -1)
    m_tabs->tabBar()->setAccessibleName(QStringLiteral("Builder pages"));
    m_tabs->addTab(helixPage(), QStringLiteral("Ideal helix"));
    m_tabs->addTab(mutatePage(), QStringLiteral("Mutate"));
    m_tabs->addTab(fitPage(), QStringLiteral("Fit to map"));
    splitter()->addWidget(m_tabs);
    QVBoxLayout* right = addReport();
    m_progress = new QProgressBar(report()->parentWidget());
    m_progress->setObjectName(nameOf(QStringLiteral("progress")));
    m_progress->setTextVisible(true);
    m_progress->setRange(0, kProgressSteps);
    m_progress->setVisible(false);
    right->addWidget(m_progress);

    m_cancel = addButton(QStringLiteral("cancel"), QStringLiteral("Cancel fit"), QStringLiteral("Stop the running density fit"));
    m_cancel->setVisible(false);
    m_check = addButton(QStringLiteral("check"), QStringLiteral("Check"),
                        QStringLiteral("Validate the inputs of the current page (nothing is written)"));
    m_write = addButton(QStringLiteral("write"), QStringLiteral("Write PDB"),
                        QStringLiteral("Run the current page and write its output file"));
    m_write->setDefault(true);  // GTK's suggested-action
    QPushButton* close = addButton(QStringLiteral("close"), QStringLiteral("Close"));
    connect(m_cancel, &QPushButton::clicked, this, &BuilderWindow::cancel);
    connect(m_check, &QPushButton::clicked, this, &BuilderWindow::checkCurrent);
    connect(m_write, &QPushButton::clicked, this, &BuilderWindow::runCurrent);
    connect(close, &QPushButton::clicked, this, &QWidget::close);
    m_jobs->addBusy(m_write);
    connect(m_tabs, &QTabWidget::currentChanged, this, &BuilderWindow::sync);
    connect(m_jobs, &jobs::JobRunner::statusChanged, this, [this](double fraction, const QString& text) {
        if (m_progress->maximum() > 0) {
            m_progress->setValue(static_cast<int>(fraction * kProgressSteps + 0.5));
            m_progress->setFormat(text);
        }
    });
    connect(m_jobs, &jobs::JobRunner::finished, this, [this](jobs::JobOutcome outcome) {
        finishJob();
        if (outcome == jobs::JobOutcome::Cancelled) showReport(m_cancelledText);
    });
    finishLayout(980, 680);
    sync();
}

BuilderWindow::~BuilderWindow() = default;

// ---- pages -----------------------------------------------------------------------------------

QWidget* BuilderWindow::helixPage() {
    QVBoxLayout* column = nullptr;
    QScrollArea* scroll = makeFormScroll(QStringLiteral("helix"), QStringLiteral("Ideal helix"), &column, m_tabs);
    QWidget* body = scroll->widget();
    const QString page = toolId() + QStringLiteral(".helix");

    auto* sequence = new pages::Section(page, QStringLiteral("sequence"), QStringLiteral("Sequence"), body);
    m_seq = makeLine(QStringLiteral("seq"), QStringLiteral("Sequence (5'->3')"), QStringLiteral("e.g. GGCAUCGAUGCC"),
                     QStringLiteral("Strand 1, 5'->3' (T is read as U)"));
    sequence->form()->add(QStringLiteral("Sequence (5'->3')"), m_seq);
    m_partnerMode = makeCombo(QStringLiteral("partner_mode"), QStringLiteral("Partner strand"),
                              {{"auto", "complement (canonical partners)"},
                               {"none", "none (single strand)"},
                               {"explicit", "explicit sequence"}},
                              "auto");
    sequence->form()->add(QStringLiteral("Partner strand"), m_partnerMode,
                          QStringLiteral("An explicit partner (5'->3', same length) may contain mismatches and G·U "
                                         "wobbles"));
    m_partner = makeLine(QStringLiteral("partner"), QStringLiteral("Partner sequence"), QStringLiteral("partner 5'->3'"));
    sequence->form()->add(QStringLiteral("Partner sequence"), m_partner);
    m_chirality = makeCombo(QStringLiteral("chirality"), QStringLiteral("Chirality"),
                            {{"D", "D-RNA (right-handed)"}, {"L", "L-RNA (mirror image, left-handed)"}}, "D");
    sequence->form()->add(QStringLiteral("Chirality"), m_chirality);
    m_pairing = makeCombo(QStringLiteral("pairing"), QStringLiteral("Pairing"),
                          {{"wc", "Watson-Crick (cWW), antiparallel A-form"},
                           {"reverse-wc", "reverse Watson-Crick (tWW), parallel"},
                           {"hoogsteen", "Hoogsteen (cWH: A·U, G·C+), parallel"},
                           {"reverse-hoogsteen", "reverse Hoogsteen (tWH), antiparallel"},
                           {"hh", "Hoogsteen/Hoogsteen A·A (tHH, PDB 4JRD), parallel"}},
                          "wc");
    sequence->form()->add(QStringLiteral("Pairing"), m_pairing,
                          QStringLiteral("Base-pair geometry of the duplex; for the Hoogsteen pairings strand 1 is the "
                                         "purine strand. Mismatches are allowed and reported."));
    m_triplex = makeCombo(QStringLiteral("triplex"), QStringLiteral("Triplex"),
                          {{"none", "none (duplex)"},
                           {"hoogsteen", "pyrimidine motif: Hoogsteen, parallel"},
                           {"reverse-hoogsteen", "purine motif: reverse Hoogsteen, antiparallel"}},
                          "none");
    sequence->form()->add(QStringLiteral("Triplex"), m_triplex,
                          QStringLiteral("A third strand in the major groove on the Hoogsteen edge of the strand-1 "
                                         "purines (Watson-Crick duplex only)"));
    m_third = makeLine(QStringLiteral("third"), QStringLiteral("Third strand"),
                       QStringLiteral("third strand 5'->3' (empty: motif bases)"),
                       QStringLiteral("Empty: U / C+ (pyrimidine motif) or A / G (purine motif) on the strand-1 A / G"));
    sequence->form()->add(QStringLiteral("Third strand"), m_third);
    column->addWidget(sequence);

    auto* atoms = new pages::Section(page, QStringLiteral("names"), QStringLiteral("Names and atoms"), body);
    m_chain1 = makeLine(QStringLiteral("chain1"), QStringLiteral("Chain of strand 1"), QStringLiteral("A"), QString(), 3);
    m_chain2 = makeLine(QStringLiteral("chain2"), QStringLiteral("Chain of strand 2"), QStringLiteral("B"), QString(), 3);
    m_chain3 = makeLine(QStringLiteral("chain3"), QStringLiteral("Chain of strand 3"), QStringLiteral("C"), QString(), 3);
    m_chain1->setText(QStringLiteral("A"));
    m_chain2->setText(QStringLiteral("B"));
    m_chain3->setText(QStringLiteral("C"));
    atoms->form()->addRow(QStringLiteral("Chains (strand 1 / 2 / 3)"), {m_chain1, m_chain2, m_chain3});
    m_start1 = makeSpin(QStringLiteral("start1"), QStringLiteral("First residue of strand 1"), -999, 9999, 1);
    m_start2 = makeSpin(QStringLiteral("start2"), QStringLiteral("First residue of strand 2"), -999, 9999, 1);
    m_start3 = makeSpin(QStringLiteral("start3"), QStringLiteral("First residue of strand 3"), -999, 9999, 1);
    atoms->form()->addRow(QStringLiteral("First residue (1 / 2 / 3)"), {m_start1, m_start2, m_start3});
    m_scheme = makeCombo(QStringLiteral("scheme"), QStringLiteral("Residue names"),
                         {{"", "default (A/C/G/U; L: 0A/0C/0G/0U)"},
                          {"amber", "Amber (L: LA/LC/LG/LU + 5'/3', Dudek & Trylska)"},
                          {"ccd", "CCD for Coot / Phenix (L: 0A … + 0A5/0A3)"},
                          {"rosetta", "Rosetta (L: 0A/0C/0G/0U, all ATOM)"},
                          {"pdb", "PDB deposition (L: 0A …, non-standard as HETATM)"}},
                         "");
    atoms->form()->add(QStringLiteral("Residue names"), m_scheme);
    m_terminals = makeCheck(QStringLiteral("terminals"), QStringLiteral("5'/3' terminal names (Amber, CCD)"), true);
    atoms->form()->addWide(m_terminals);
    m_phosphate5 = makeCheck(QStringLiteral("phosphate5"), QStringLiteral("Keep a 5'-phosphate (P, OP1, OP2, OP3)"));
    atoms->form()->addWide(m_phosphate5);
    m_hydrogens = makeCheck(QStringLiteral("hydrogens"), QStringLiteral("Add hydrogens (default: heavy atoms only)"));
    atoms->form()->addWide(m_hydrogens);
    column->addWidget(atoms);

    auto* output = new pages::Section(page, QStringLiteral("output"), QStringLiteral("Output"), body);
    m_helixOut = makePath(QStringLiteral("helix_out"), QStringLiteral("Output PDB"), QStringLiteral("Output PDB"),
                          forms::PathMode::SaveFile, pdbFiles());
    output->form()->add(QStringLiteral("Output PDB"), m_helixOut);
    column->addWidget(output);
    column->addWidget(pages::makeHint(
        QStringLiteral("Watson-Crick: Arnott A-RNA fiber geometry, C3'-endo sugars, 11 bp per turn (twist 32.7°, rise "
                       "2.81 Å). Other pairings and triplexes: the PDB 4JRD nucleotide with fitted helices and backbones "
                       "(starting models: minimise; the report lists link gaps, clashes and mismatches). The L helix is "
                       "the exact mirror image of the D helix."),
        body));
    column->addStretch(1);

    for (QComboBox* combo : {m_partnerMode, m_pairing, m_triplex}) {
        connect(combo, &QComboBox::currentIndexChanged, this, &BuilderWindow::sync);
    }
    return scroll;
}

QWidget* BuilderWindow::mutatePage() {
    QVBoxLayout* column = nullptr;
    QScrollArea* scroll = makeFormScroll(QStringLiteral("mutate"), QStringLiteral("Mutate"), &column, m_tabs);
    QWidget* body = scroll->widget();
    auto* box = new pages::Section(toolId() + QStringLiteral(".mutate"), QStringLiteral("mutations"),
                                   QStringLiteral("Mutations"), body);
    m_mutInput = makePath(QStringLiteral("mut_input"), QStringLiteral("Input PDB"), QStringLiteral("Select a PDB file"),
                          forms::PathMode::OpenFile, pdbFiles());
    box->form()->add(QStringLiteral("Input PDB"), m_mutInput);
    m_mutations = makeLine(QStringLiteral("mutations"), QStringLiteral("Mutations"), QStringLiteral("e.g. A:5G, B:20C"));
    box->form()->add(QStringLiteral("Mutations"), m_mutations, QStringLiteral("Chain:residue number followed by the new base"));
    m_compensatory = makeCheck(QStringLiteral("compensatory"), QStringLiteral("Compensatory mutations of the paired partners"),
                               true,
                               QStringLiteral("Every mutated base that is Watson-Crick / G·U paired in the input gets its "
                                              "partner changed to the complement, unless you mutate the partner yourself"));
    box->form()->addWide(m_compensatory);
    m_wobble = makeCheck(QStringLiteral("wobble"), QStringLiteral("Keep G·U pairs where possible"));
    box->form()->addWide(m_wobble);
    m_mutOut = makePath(QStringLiteral("mut_out"), QStringLiteral("Output PDB"), QStringLiteral("Output PDB"),
                        forms::PathMode::SaveFile, pdbFiles());
    box->form()->add(QStringLiteral("Output PDB"), m_mutOut);
    column->addWidget(box);
    column->addWidget(pages::makeHint(
        QStringLiteral("The new base is superposed on the old glycosidic frame (C1', N9/N1 and its two ring neighbours); "
                       "sugar and backbone atoms are kept as they are. Residue names keep their scheme and chirality "
                       "(LG -> LC, 0G -> 0C)."),
        body));
    column->addStretch(1);
    connect(m_mutInput, &PathEntry::pathChanged, this, [this] { suggest(QStringLiteral("mutate")); });
    connect(m_compensatory, &QCheckBox::toggled, this, &BuilderWindow::sync);
    return scroll;
}

QWidget* BuilderWindow::fitPage() {
    QVBoxLayout* column = nullptr;
    QScrollArea* scroll = makeFormScroll(QStringLiteral("fit"), QStringLiteral("Fit to map"), &column, m_tabs);
    QWidget* body = scroll->widget();
    auto* box = new pages::Section(toolId() + QStringLiteral(".fit"), QStringLiteral("fit"),
                                   QStringLiteral("Rigid-body fit"), body);
    m_fitInput = makePath(QStringLiteral("fit_input"), QStringLiteral("Model PDB"), QStringLiteral("Select a PDB file"),
                          forms::PathMode::OpenFile, pdbFiles());
    box->form()->add(QStringLiteral("Model PDB"), m_fitInput);
    m_fitMap = makePath(QStringLiteral("fit_map"), QStringLiteral("Density map"), QStringLiteral("Select a density map"),
                        forms::PathMode::OpenFile, mapFiles());
    box->form()->add(QStringLiteral("Density map"), m_fitMap);
    m_resolution = makeDoubleSpin(QStringLiteral("resolution"), QStringLiteral("Resolution (Å)"), 0.5, 50.0, 0.1, 1, 4.0);
    box->form()->add(QStringLiteral("Resolution (Å)"), m_resolution, QStringLiteral("Width of the map simulated from the model"));
    m_search = makeCombo(QStringLiteral("search"), QStringLiteral("Search"),
                         {{"global", "global search (all orientations)"},
                          {"local", "local refinement of the model as placed"}},
                         "global");
    box->form()->add(QStringLiteral("Search"), m_search);
    m_step = makeDoubleSpin(QStringLiteral("step"), QStringLiteral("Angular step (°)"), 5.0, 45.0, 5.0, 0, 30.0);
    box->form()->add(QStringLiteral("Angular step (°)"), m_step);
    m_fitOut = makePath(QStringLiteral("fit_out"), QStringLiteral("Output PDB"), QStringLiteral("Output PDB"),
                        forms::PathMode::SaveFile, pdbFiles());
    box->form()->add(QStringLiteral("Output PDB"), m_fitOut);
    column->addWidget(box);
    column->addWidget(
        pages::makeHint(QStringLiteral("Only coordinates change: every other line of the model is written as it was."), body));
    column->addStretch(1);
    connect(m_fitInput, &PathEntry::pathChanged, this, [this] { suggest(QStringLiteral("fit")); });
    connect(m_search, &QComboBox::currentIndexChanged, this, &BuilderWindow::sync);
    return scroll;
}

void BuilderWindow::sync() {
    const QString mode = currentId(m_partnerMode);
    m_partner->setEnabled(mode == QLatin1String("explicit"));
    const bool duplex = mode != QLatin1String("none");
    m_chain2->setEnabled(duplex);
    m_start2->setEnabled(duplex);
    m_pairing->setEnabled(duplex);
    m_triplex->setEnabled(duplex && currentId(m_pairing) == QLatin1String("wc"));
    const bool triple = m_triplex->isEnabled() && currentId(m_triplex) != QLatin1String("none");
    for (QWidget* w : std::initializer_list<QWidget*>{m_third, m_chain3, m_start3}) w->setEnabled(triple);
    m_wobble->setEnabled(m_compensatory->isChecked());
    m_step->setEnabled(currentId(m_search) == QLatin1String("global"));
    m_check->setEnabled(page() != QLatin1String("fit"));
}

// <input>_mut.pdb / <input>_fit.pdb while the output is empty or still the previous suggestion
void BuilderWindow::suggest(const QString& which) {
    const bool mut = which == QLatin1String("mutate");
    PathEntry* src = mut ? m_mutInput : m_fitInput;
    PathEntry* dst = mut ? m_mutOut : m_fitOut;
    QString& automatic = mut ? m_autoMutate : m_autoFit;
    const QString current = dst->path();
    if (!current.isEmpty() && current != automatic) {
        return;
    }
    automatic = toQString(suggestedOutput(toStdString(src->path()), mut ? "mut" : "fit"));
    dst->setPath(automatic);
}

// ---- actions ---------------------------------------------------------------------------------

QString BuilderWindow::page() const {
    const int k = m_tabs->currentIndex();
    return QString::fromLatin1(kPages[k >= 0 && k < 3 ? k : 0]);
}

void BuilderWindow::setPage(const QString& name) {
    for (int k = 0; k < 3; ++k) {
        if (name == QLatin1String(kPages[k])) m_tabs->setCurrentIndex(k);
    }
}

bool BuilderWindow::busy() const { return m_jobs->running(); }

HelixRequest BuilderWindow::helixRequest() const {
    HelixRequest r;
    const QString mode = currentId(m_partnerMode);
    r.partner = mode == QLatin1String("explicit") ? toStdString(m_partner->text().trimmed()) : toStdString(mode);
    const QString c1 = m_chain1->text().trimmed();
    const QString c2 = m_chain2->text().trimmed();
    r.pairing = mode != QLatin1String("none") ? toStdString(currentId(m_pairing)) : std::string("wc");
    const std::string triplex =
        mode != QLatin1String("none") && r.pairing == "wc" ? toStdString(currentId(m_triplex)) : std::string("none");
    r.chains = {c1.isEmpty() ? std::string("A") : toStdString(c1), c2.isEmpty() ? std::string("B") : toStdString(c2)};
    r.start = {m_start1->value(), m_start2->value()};
    if (triplex != "none") {
        const QString c3 = m_chain3->text().trimmed();
        r.chains.push_back(c3.isEmpty() ? std::string("C") : toStdString(c3));
        r.start.push_back(m_start3->value());
    }
    r.sequence = toStdString(m_seq->text());
    r.chirality = toStdString(currentId(m_chirality));
    r.scheme = toStdString(currentId(m_scheme));
    r.terminals = m_terminals->isChecked();
    r.hydrogens = m_hydrogens->isChecked();
    r.phosphate5 = m_phosphate5->isChecked();
    r.triplex = triplex == "none" ? std::string() : triplex;
    const QString third = m_third->text().trimmed();
    r.third = third.isEmpty() ? std::string("auto") : toStdString(third);
    r.output = toStdString(m_helixOut->path());
    return r;
}

QString BuilderWindow::inputFile(PathEntry* entry, const QString& what) {
    const QString path = entry->path();
    if (path.isEmpty() || !isFile(toStdString(path))) {
        showReport(QStringLiteral("Choose an existing %1 first.").arg(what));
        return {};
    }
    return path;
}

QString BuilderWindow::outputFile(PathEntry* entry, const QString& source) {
    const QString out = entry->path();
    if (out.isEmpty()) {
        showReport(QStringLiteral("Choose an output file."));
        return {};
    }
    if (!source.isEmpty() && samePath(toStdString(out), toStdString(source))) {
        showReport(QStringLiteral("The output must not overwrite the input file."));
        return {};
    }
    return out;
}

void BuilderWindow::checkCurrent() {
    const QString current = page();
    if (current == QLatin1String("helix")) {
        try {
            showReport(toQString(helixCheckReport(helixRequest())));
        } catch (const std::exception& e) {
            showReport(QStringLiteral("ERROR: ") + toQString(errorText(e)));
        }
    } else if (current == QLatin1String("mutate")) {
        const QString path = inputFile(m_mutInput, QStringLiteral("input PDB"));
        if (path.isEmpty()) {
            return;
        }
        try {
            showReport(toQString(mutateCheckReport(toStdString(path), toStdString(m_mutations->text()))));
        } catch (const std::exception& e) {  // TW-02 (MDATK 0.1.0): also a file that is not UTF-8 text
            showReport(QStringLiteral("ERROR: ") + toQString(errorText(e)));
        }
    }
}

void BuilderWindow::runCurrent() {
    if (busy()) {
        return;  // MDATK 0.1.1: a job of the worker has not reported back yet
    }
    const QString current = page();
    if (current == QLatin1String("helix")) {
        runHelix();
    } else if (current == QLatin1String("mutate")) {
        runMutate();
    } else {
        runFit();
    }
}

void BuilderWindow::cancel() { m_jobs->cancel(); }

bool BuilderWindow::mayClose() { return m_jobs->confirmClose(this); }

void BuilderWindow::startJob(const QString& message, bool cancellable, bool pulse) {
    showReport(message);
    m_progress->setRange(0, pulse ? 0 : kProgressSteps);  // 0, 0: busy (MDATK pulses the bar)
    m_progress->setValue(0);
    m_progress->setFormat(QString());
    m_progress->setVisible(true);
    m_cancel->setVisible(cancellable);
}

void BuilderWindow::finishJob() {
    m_progress->setVisible(false);
    m_progress->setRange(0, kProgressSteps);
    m_cancel->setVisible(false);
}

void BuilderWindow::showWritten(const ToolResult& result) {
    QString head = QStringLiteral("Wrote ") + gui::displayPath(toQString(result.output));
    if (result.cc) {
        char cc[32];
        std::snprintf(cc, sizeof cc, "%.3f", *result.cc);
        head += QStringLiteral(" (CC %1)").arg(QString::fromLatin1(cc));
    }
    QStringList lines{head, QString()};
    for (const std::string& note : result.notes) lines << toQString(note);
    showReport(lines.join(QLatin1Char('\n')));
    QString log = QStringLiteral("Structure builder: wrote ") + gui::displayPath(toQString(result.output));
    if (result.cc) {
        char cc[32];
        std::snprintf(cc, sizeof cc, "%.3f", *result.cc);
        log += QStringLiteral(" (density fit CC %1)").arg(QString::fromLatin1(cc));
    }
    host().appendLog(log);
}

void BuilderWindow::runHelix() {
    const QString out = outputFile(m_helixOut);
    if (out.isEmpty()) {
        return;
    }
    HelixRequest r = helixRequest();
    r.output = toStdString(out);
    const fs::path project = host().projectPath();
    if (!helixInBackground(r)) {
        try {
            showWritten(buildHelix(r, project));
        } catch (const std::exception& e) {
            showReport(QStringLiteral("ERROR: ") + toQString(errorText(e)));
        }
        return;
    }
    // the fitted pairings can take tens of seconds: build in the background
    startJob(QStringLiteral("Building the helix (other pairings take a few seconds per 10 nucleotides) …"), false, true);
    m_cancelledText = QStringLiteral("Cancelled: the helix build stopped, nothing was written.");
    m_jobs->start(
        QStringLiteral("Building the helix"),
        [r, project](const core::Reporter& reporter) { return buildHelix(r, project, reporter); },
        [this](const ToolResult& result) { showWritten(result); },
        [this](const jobs::JobFailure& failure) { showReport(QStringLiteral("ERROR: ") + failureText(failure)); });
}

void BuilderWindow::runMutate() {
    const QString path = inputFile(m_mutInput, QStringLiteral("input PDB"));
    if (path.isEmpty()) {
        return;
    }
    const QString out = outputFile(m_mutOut, path);
    if (out.isEmpty()) {
        return;
    }
    try {
        showWritten(mutate(toStdString(path), toStdString(m_mutations->text()), m_compensatory->isChecked(),
                           m_wobble->isChecked(), toStdString(out), host().projectPath()));
    } catch (const std::exception& e) {
        showReport(QStringLiteral("ERROR: ") + toQString(errorText(e)));
    }
}

void BuilderWindow::runFit() {
    if (busy()) {
        return;
    }
    const QString path = inputFile(m_fitInput, QStringLiteral("model PDB"));
    if (path.isEmpty()) {
        return;
    }
    const QString map = inputFile(m_fitMap, QStringLiteral("density map"));
    if (map.isEmpty()) {
        return;
    }
    const QString out = outputFile(m_fitOut, path);
    if (out.isEmpty()) {
        return;
    }
    const std::string input = toStdString(path);
    const std::string mapPath = toStdString(map);
    const std::string output = toStdString(out);
    const double resolution = m_resolution->value();
    const std::string search = toStdString(currentId(m_search));
    const double step = m_step->value();
    const fs::path project = host().projectPath();
    startJob(QStringLiteral("Fitting %1 into %2 …").arg(basename(path), basename(map)), true, false);
    // TW-03 (MDATK 0.0.13): a cancelled fit writes nothing and is logged as cancelled
    m_cancelledText = QStringLiteral("Cancelled: the density fit stopped, nothing was written.");
    m_jobs->start(
        QStringLiteral("Fitting %1 into %2").arg(basename(path), basename(map)),
        [=](const core::Reporter& reporter) {
            return fitToMap(input, mapPath, resolution, search, step, output, project, reporter);
        },
        [this](const ToolResult& result) { showWritten(result); },
        [this](const jobs::JobFailure& failure) { showReport(QStringLiteral("ERROR: ") + failureText(failure)); });
}

}  // namespace mdbridge::app::tools
