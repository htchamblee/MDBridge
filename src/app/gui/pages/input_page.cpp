#include "gui/pages/input_page.hpp"

#include <string>
#include <vector>

#include <QHBoxLayout>
#include <QPlainTextEdit>
#include <QPushButton>
#include <QLabel>
#include <QTableView>

#include "gui/forms/field_widgets.hpp"
#include "gui/forms/list_editor.hpp"
#include "mdbridge/core/conda.hpp"

namespace mdbridge::app::pages {

namespace {

void setTextHeight(forms::FieldWidget* field, int height) {
    if (auto* text = qobject_cast<QPlainTextEdit*>(field->editor())) {
        text->setMinimumHeight(height);
        text->setMaximumHeight(height * 2);
    }
}

void setTableHeight(forms::FieldWidget* field, int height) {
    if (auto* list = qobject_cast<forms::ListEditor*>(field)) {
        list->view()->setMinimumHeight(height);
    }
}

}  // namespace

InputPage::InputPage(forms::FormContext& context, QWidget* parent)
    : FormPage(context, QStringLiteral("input"), QStringLiteral("Input & cleaning"), parent) {
    // ---- Project notes
    Section* notes = addSection(QStringLiteral("notes"), tr("Project notes"));
    forms::FieldWidget* text = fieldAs("notes", tr("Project notes"));
    setTextHeight(text, 54);
    notes->form()->addWide(text);
    notes->form()->addWide(makeHint(tr("Free text saved in the project file (aims, decisions, to-dos). The steps you "
                                "run are logged there too (menu → Project log).")));

    // ---- Input files
    Section* files = addSection(QStringLiteral("files"), tr("Input files"));
    files->form()->add(tr("Topology"), fieldAs("inputs.topology", tr("Topology")),
               tr("Amber .parm7 (preferred), cleaned PDB or another topology"));
    forms::FieldWidget* trajectories = fieldAs("inputs.trajectories", tr("Trajectories"));
    setTableHeight(trajectories, 90);
    files->form()->add(tr("Trajectories"), trajectories,
               tr("Amber NetCDF .nc (preferred), DCD, XTC, ... (several files are concatenated)"));

    // ---- Frames, time and external data
    Section* frames = addSection(QStringLiteral("frames"), tr("Frames, time and external data"));
    const QString frameHint = tr("0-based frames of the (cleaned) trajectory; stop is exclusive (blank = last)");
    frames->form()->addRow(tr("Start / stop / step"),
                   {fieldAs("inputs.start", tr("Start"), frameHint), fieldAs("inputs.stop", tr("Stop"), frameHint),
                    fieldAs("inputs.step", tr("Step"), frameHint)},
                   frameHint);
    frames->form()->add(tr("Time step override (ps)"), fieldAs("inputs.dt_ps", tr("Time step override (ps)")),
                tr("Leave blank to use the trajectory's time stamps"));
    frames->form()->add(tr("Display time unit"), fieldAs("inputs.time_unit", tr("Display time unit")));
    forms::FieldWidget* external = fieldAs("inputs.external", tr("External data (CSV/text)"));
    setTableHeight(external, 90);
    frames->form()->add(tr("External data (CSV/text)"), external,
                tr("Extra traces (e.g. cpptraj .dat files) interpolated onto the analysed frames"));

    // ---- Trajectory cleaning (cpptraj)
    Section* cl = addSection(QStringLiteral("cleaning"), tr("Trajectory cleaning (cpptraj)"));
    cl->form()->addWide(fieldAs("cleaning.enabled", tr("Strip and clean the trajectory with cpptraj before analysis")));
    cl->form()->add(tr("cpptraj executable"), fieldAs("cleaning.cpptraj_path", tr("cpptraj executable")),
            tr("Name on PATH or absolute path"));
    m_detect = makeButton(tr("Detect"), QStringLiteral("cleaning.detect"),
                          tr("Find environments that contain bin/cpptraj"));
    cl->form()->addRow(tr("Conda environment"), {fieldAs("cleaning.conda_env", tr("Conda environment")), m_detect},
               tr("Activated for every cpptraj call (conda run -n <env> cpptraj); a name or a path"));
    cl->form()->add(tr("conda executable"), fieldAs("cleaning.conda_exe", tr("conda executable")),
            tr("Blank: found automatically (conda, mamba or micromamba)"));
    cl->form()->add(tr("Solvent mask"), fieldAs("cleaning.solvent_mask", tr("Solvent mask")));
    cl->form()->add(tr("Ion mask"), fieldAs("cleaning.ion_mask", tr("Ion mask")), tr("Mg2+ is kept by default"));
    cl->form()->add(tr("Extra strip mask"), fieldAs("cleaning.extra_strip_mask", tr("Extra strip mask")));
    cl->form()->addRow(tr("Imaging"), {fieldAs("cleaning.autoimage", tr("autoimage, anchor")),
                               fieldAs("cleaning.autoimage_anchor", tr("Autoimage anchor"))});
    cl->form()->addRow(tr("Fitting"), {fieldAs("cleaning.fit", tr("rms fit, mask")), fieldAs("cleaning.fit_mask", tr("Fit mask"))});
    cl->form()->addRow(tr("Fit reference"), {fieldAs("cleaning.fit_reference", tr("Fit reference")),
                                     fieldAs("cleaning.reference_frame", tr("Reference frame"))});
    cl->form()->add(tr("Reference file"), fieldAs("cleaning.reference_file", tr("Reference file")));
    cl->form()->addWide(fieldAs("cleaning.remove_box_after_fit",
                      tr("Remove box information after fitting (a rotated box cannot be stored in NetCDF)")));
    const QString trajinHint = tr("1-based cpptraj trajin arguments");
    cl->form()->addRow(tr("trajin start / stop / offset"),
               {fieldAs("cleaning.start", tr("trajin start"), trajinHint), fieldAs("cleaning.stop", tr("trajin stop"), trajinHint),
                fieldAs("cleaning.offset", tr("trajin offset"), trajinHint)},
               trajinHint);
    forms::FieldWidget* extra = fieldAs("cleaning.extra_commands", tr("Extra cpptraj commands"));
    setTextHeight(extra, 50);
    cl->form()->add(tr("Extra cpptraj commands"), extra, tr("Advanced: inserted before trajout"));
    auto* buttons = new QWidget;
    auto* row = new QHBoxLayout(buttons);
    row->setContentsMargins(0, 0, 0, 0);
    m_runCleaning = makeButton(tr("Run cleaning now"), QStringLiteral("cleaning.run"));
    m_runCleaning->setIcon(style()->standardIcon(QStyle::SP_MediaPlay));
    m_inspect = makeButton(tr("Inspect system"), QStringLiteral("cleaning.inspect"),
                           tr("Print atom/residue/frame counts and check residue groups"));
    row->addWidget(m_runCleaning);
    row->addWidget(m_inspect);
    row->addStretch(1);
    cl->form()->addWide(buttons);

    connect(m_detect, &QPushButton::clicked, this, [this] { detectConda(); });
    connect(m_runCleaning, &QPushButton::clicked, this, &InputPage::runCleaningRequested);
    connect(m_inspect, &QPushButton::clicked, this, &InputPage::inspectSystemRequested);
}

QString InputPage::detectConda() {
    const std::vector<core::conda::CondaEnv> envs = core::conda::findCpptrajEnvs();
    if (envs.empty()) {
        std::string roots;
        for (const std::string& r : core::conda::candidateRoots()) {
            roots += (roots.empty() ? "" : ", ") + r;
        }
        emit information(tr("No conda environment with cpptraj found"),
                         tr("Searched installations: %1. Enter the environment name or path manually.")
                             .arg(roots.empty() ? tr("none found") : forms::toQString(roots)));
        return {};
    }
    const QString name = forms::toQString(envs.front().name);
    context().document().setText("cleaning.conda_env", envs.front().name, "Detect conda environment");
    std::string found;
    for (const core::conda::CondaEnv& e : envs) {
        found += (found.empty() ? "" : ", ") + e.name + " (" + e.prefix + ")";
    }
    emit logMessage(tr("Conda environments with cpptraj: %1").arg(forms::toQString(found)));
    return name;
}

forms::FieldWidget* InputPage::fieldAs(const std::string& path, const QString& label, const QString& hint) {
    return field(path, [&](forms::FieldSpec& spec) {
        spec.label = label;
        if (spec.hint.isEmpty()) {
            spec.hint = hint;
        }
    });
}

}  // namespace mdbridge::app::pages
