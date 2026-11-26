#include "gui/pages/outputs_page.hpp"

#include <optional>
#include <string>
#include <string_view>
#include <vector>

#include <QLabel>
#include <QPushButton>

#include "gui/forms/custom_fields.hpp"
#include "mdbridge/render/panels/charts.hpp"

namespace mdbridge::app::pages {

namespace {

QLabel* times() { return new QLabel(QStringLiteral("×")); }

// The formats offered as boxes: every format the charts write (RN-11), without the aliases jpeg and
// tiff (a list that has them gets a box of its own).
forms::FieldSpec formatSpec() {
    forms::FieldSpec spec;
    spec.label = QStringLiteral("Formats");
    for (std::string_view f : render::charts::chartFormats()) {
        if (f == "jpeg" || f == "tiff") {
            continue;
        }
        spec.choices.push_back({forms::toQString(f), {}});
    }
    spec.hint = QStringLiteral("Files written for each chart; none ticked: PNG");
    return spec;
}

}  // namespace

OutputsPage::OutputsPage(forms::FormContext& context, QWidget* parent)
    : FormPage(context, QStringLiteral("outputs"), tr("Outputs"), parent) {
    buildLocations();
    buildVideo();
    buildCharts();
    buildPdbs();
}

void OutputsPage::buildLocations() {
    Section* s = addSection(QStringLiteral("locations"), tr("Locations"));
    FormGrid* f = s->form();
    f->add(tr("Output directory"), field("output.out_dir"));
    f->add(tr("Cache directory"), field("output.cache_dir"));
    // MDATK's check box text wraps; a Qt check box does not, so its second half is a hint below it
    f->addWide(check("output.write_params_log", tr("Write MDATK_params_log.toml"),
                     tr("Parameters, versions, cache keys and citations, for reproducibility")));
    f->addHint(tr("(parameters, versions, cache keys, citations) for reproducibility; open it as a project to "
                  "repeat the run"));
}

void OutputsPage::buildVideo() {
    Section* s = addSection(QStringLiteral("video"), tr("Video"));
    FormGrid* f = s->form();
    f->addWide(check("output.video.enabled", tr("Write MP4 video")));
    f->add(tr("File name"), field("output.video.filename"));
    f->addRow(tr("Resolution (px)"), {field("output.video.width"), times(), field("output.video.height")});
    f->add(tr("DPI"), field("output.video.dpi"));
    f->addRow(tr("FPS / frame stride"), {field("output.video.fps"), field("output.video.frame_stride")});
    f->addRow(tr("Codec / quality (CRF)"), {field("output.video.codec"), field("output.video.crf")});
    f->add(tr("FFmpeg executable"), field("output.video.ffmpeg_path"));
}

void OutputsPage::buildCharts() {
    Section* s = addSection(QStringLiteral("charts"), tr("Time charts"));
    FormGrid* f = s->form();
    f->addWide(check("output.charts.enabled", tr("Write standalone time charts")));
    f->add(tr("Series"), field("output.charts.series"));
    f->addWide(check("output.charts.include_heatmaps", tr("Also chart per-residue heatmaps")));
    m_formats = new forms::ChoiceChecksField(
        context(), "output.charts.formats", formatSpec(),
        [](const QString& v) { return forms::toQString(render::charts::normalizeFormat(forms::toStdString(v))); },
        [](const QString& v) { return !render::charts::chartFormatProblem({forms::toStdString(v)}).has_value(); },
        [](const QStringList& list) {
            std::vector<std::string> formats;
            for (const QString& v : list) {
                formats.push_back(forms::toStdString(v));
            }
            const std::optional<std::string> problem = render::charts::chartFormatProblem(formats);
            return problem ? forms::toQString(*problem) : QString();
        });
    f->add(tr("Formats"), m_formats);
    f->addRow(tr("Size (in)"), {field("output.charts.width"), times(), field("output.charts.height")});
    f->add(tr("DPI"), field("output.charts.dpi"));
    f->addRow(tr("Style / line colour / width"),
              {field("output.charts.style",
                     [](forms::FieldSpec& spec) {
                         // AP-19: the styles MDBridge draws (owner decision D5 of milestone 4); a typed
                         // name is kept and warned about when the charts are drawn
                         spec.choices = {{QStringLiteral("mdatk"), {}}, {QStringLiteral("default"), {}}};
                         spec.hint = QStringLiteral("mdatk: MDATK's chart look; default: Matplotlib's. Other "
                                                    "names are kept and drawn with MDATK's look.");
                     }),
               field("output.charts.color"), field("output.charts.line_width")});
    f->add(tr("Title template"), field("output.charts.title"));
    f->addRow(tr("Axis labels"), {field("output.charts.xlabel"), field("output.charts.ylabel")});
    f->addRow(tr("Smoothing (frames)"),
              {field("output.charts.smoothing"), check("output.charts.show_target", tr("show target line"))});
    f->add(tr("File prefix"), field("output.charts.prefix"));
}

void OutputsPage::buildPdbs() {
    Section* s = addSection(QStringLiteral("pdb"), tr("Frame PDBs"));
    FormGrid* f = s->form();
    f->addWide(check("output.pdb.enabled", tr("Export frames as PDB")));
    f->add(tr("Frames"), field("output.pdb.frames"));
    f->add(tr("Times"), field("output.pdb.times", [](forms::FieldSpec& spec) {
               spec.hint = QStringLiteral("Nearest frame is used; unitless values use the display time unit");
           }));
    f->add(tr("Atom subset"), field("output.pdb.selection", [](forms::FieldSpec& spec) {
               spec.hint = QStringLiteral("Group names (comma separated) or an atom selection");
           }));
    f->add(tr("File prefix"), field("output.pdb.prefix"));
    m_exportPdbs = makeButton(tr("Export PDBs now"), QStringLiteral("output.export_pdbs"));
    m_rerender = makeButton(tr("Re-render video + charts from cache"), QStringLiteral("output.rerender"));
    f->addButtons({m_exportPdbs, m_rerender});
    connect(m_exportPdbs, &QPushButton::clicked, this, [this] {
        context().commitPendingEdits();
        emit exportPdbsRequested();
    });
    connect(m_rerender, &QPushButton::clicked, this, [this] {
        context().commitPendingEdits();
        emit rerenderRequested();
    });
    setJobActionsAvailable(false);
}

void OutputsPage::setJobActionsAvailable(bool available) {
    m_exportPdbs->setEnabled(available);
    m_rerender->setEnabled(available);
    const QString rerender = tr("Uses the last cached analyses (no trajectory access); handy after changing "
                                "colours, layout or chart style");
    m_exportPdbs->setToolTip(available ? QString() : tr("Available in the analysis window"));
    m_rerender->setToolTip(available ? rerender : rerender + tr("\nAvailable in the analysis window"));
}

}  // namespace mdbridge::app::pages
