#include "gui/pages/structure_page.hpp"

#include <algorithm>
#include <string>

#include <QLabel>
#include <QPointer>
#include <QPushButton>
#include <QTableView>
#include <QTimer>

#include "gui/forms/custom_fields.hpp"
#include "gui/forms/list_editor.hpp"
#include "mdbridge/core/project_document.hpp"

namespace mdbridge::app::pages {

namespace {

bool hasPlacement(const core::DensityParams& d) {
    return !d.transform.empty() || !d.anchor.empty() || d.fit_cc != 0.0 || !d.fit_info.empty();
}

void clearPlacement(core::DensityParams& d) {
    d.transform.clear();
    d.anchor.clear();
    d.fit_cc = 0.0;
    d.fit_info.clear();
}

}  // namespace

StructurePage::StructurePage(forms::FormContext& context, QWidget* parent)
    : FormPage(context, QStringLiteral("structure"), tr("Structure view"), parent) {
    buildOrientation();
    buildRepresentations();
    buildDensity();
    buildTags();
    buildAppearance();
    connect(&context, &forms::FormContext::documentChanged, this,
            [this](const core::DocumentChange& change) { documentChanged(change); });
    m_densityStatus->setText(densityStatusText());
}

void StructurePage::buildOrientation() {
    Section* s = addSection(QStringLiteral("orientation"), tr("Orientation"));
    FormGrid* f = s->form();
    m_rotation = new forms::RotationField(context(), "visual.rotation", forms::defaultFieldSpec("visual.rotation"));
    f->addWide(m_rotation);
    f->add(tr("Zoom"), field("visual.zoom", [](forms::FieldSpec& spec) {
               spec.hint = QStringLiteral("1 fits the whole model; 2 shows it twice as large (mouse wheel over "
                                          "the structure panel)");
           }));
    f->add(tr("Centre on"), field("visual.center_on", [](forms::FieldSpec& spec) {
               spec.hint = QStringLiteral("Residue (number, or chain:number) or 'origin' (0, 0, 0 of the "
                                          "coordinates) in the middle of the view; empty: the model centre");
           }));
    f->addWide(check("visual.center_follow", tr("Follow the centre residue in every frame")));
    QPushButton* resetOrientation = makeButton(tr("Reset orientation"), QStringLiteral("structure.reset_orientation"));
    QPushButton* resetView = makeButton(tr("Reset zoom and centre"), QStringLiteral("structure.reset_view"));
    f->addButtons({resetOrientation, resetView});
    connect(resetOrientation, &QPushButton::clicked, this, [this] {
        context().document().set("visual.rotation", core::Value(core::Array{0.0, 0.0, 0.0}), "Reset orientation");
    });
    connect(resetView, &QPushButton::clicked, this, [this] {
        context().commitPendingEdits();
        context().document().edit("Reset zoom and centre", [](core::Project& p) {
            p.visual.zoom = 1.0;
            p.visual.center_on.clear();
        });
    });
    f->addHint(tr("Drag the structure panel in the figure to rotate the model, use the mouse wheel over it "
                  "to zoom. Orientation, zoom and centre are set on the reference frame and applied to every "
                  "frame of the video."));
}

void StructurePage::buildRepresentations() {
    Section* s = addSection(QStringLiteral("representations"), tr("Representations (3D)"));
    FormGrid* f = s->form();
    m_reps = qobject_cast<forms::ListEditor*>(field("visual.representations"));
    m_reps->addButton()->setText(tr("Add representation"));
    m_reps->addButton()->setToolTip(tr("Draw residue groups as cartoon, ribbon, tube, licorice, VDW or surface"));
    m_reps->removeButton()->setToolTip(tr("Remove the selected representation"));
    m_reps->view()->setMinimumHeight(110);
    m_reps->view()->setMaximumHeight(170);
    f->addWide(m_reps);
    QPushButton* apply = makeButton(tr("Apply"), QStringLiteral("structure.apply"), tr("Redraw the preview"));
    connect(apply, &QPushButton::clicked, this, [this] {
        context().commitPendingEdits();
        emit previewRequested();
    });
    f->addButtons({apply});
    const QString restHint = tr("Residues no representation draws are shown as a 3D trace or hidden. Ribbons, "
                                "cartoons and the trace run through P or C4'.");
    f->addRow(tr("Other residues / backbone atom"),
              {field("visual.rep_rest"), field("visual.backbone_atom"), check("visual.rep_hydrogens", tr("hydrogens"))},
              restHint);
    f->add(tr("Renderer"), field("visual.renderer", [](forms::FieldSpec& spec) {
               spec.hint = QStringLiteral("The built-in renderer needs no extra packages and works without a "
                                          "display (clusters). The GPU renderer (Vulkan) draws the same scene "
                                          "faster; where it cannot render, the built-in one is used.");
           }));
    f->addHint(tr("No representation: the 2D backbone trace. Each representation draws its groups (empty: all "
                  "shown groups) and they intersect by depth. Colouring 'auto' uses the residue colours (groups "
                  "/ proximity / glow) for backbones, nucleotide colours for cartoon bases and element colours "
                  "for atoms."));
}

void StructurePage::buildDensity() {
    Section* s = addSection(QStringLiteral("density"), tr("Density map"));
    FormGrid* f = s->form();
    f->add(tr("Map file"), field("density.map_path"));
    f->addRow(tr("Show / colour / opacity"),
              {check("density.show", tr("show")), field("density.color"), field("density.opacity"),
               check("density.front", tr("in front"),
                     tr("Draw the shadow over the 3D representations instead of behind them"))},
              tr("The map is drawn as a projected shadow behind the model"));
    f->addRow(tr("Contour level (σ) / resolution (Å)"), {field("density.level"), field("density.resolution")},
              tr("Density below mean + level × σ casts no shadow; the resolution sets the width of the "
                 "simulated map used for fitting"));
    f->addRow(tr("Fit frame / start / step (°)"),
              {field("density.fit_frame"), field("density.search"), field("density.angular_step")},
              tr("Frame 'ref' is the first analysed frame (the video's reference); any trajectory frame number "
                 "works. The global search samples orientations every <step> degrees."));
    f->addWide(check("density.auto_fit", tr("Fit automatically before rendering when the map has no placement yet")));

    m_globalSearch = makeButton(tr("Global search"), QStringLiteral("density.global_search"),
                                tr("Exhaustive orientation search + local refinement (maximises the "
                                   "correlation with a map simulated from the model)"));
    m_refine = makeButton(tr("Refine"), QStringLiteral("density.refine"),
                          tr("Local refinement from the current placement"));
    m_centre = makeButton(tr("Centre on model"), QStringLiteral("density.centre_on_model"));
    m_clear = makeButton(tr("Clear"), QStringLiteral("density.clear"), tr("Forget the placement"));
    setCentreOnModelAvailable(false);
    f->addButtons({m_globalSearch, m_refine, m_centre, m_clear});
    connect(m_globalSearch, &QPushButton::clicked, this, [this] {
        context().commitPendingEdits();
        emit fitRequested(QStringLiteral("global"));
    });
    connect(m_refine, &QPushButton::clicked, this, [this] {
        context().commitPendingEdits();
        emit fitRequested(QStringLiteral("manual"));
    });
    connect(m_centre, &QPushButton::clicked, this, [this] {
        context().commitPendingEdits();
        emit centreOnModelRequested();
    });
    connect(m_clear, &QPushButton::clicked, this, [this] {
        context().document().edit("Clear the map placement", [](core::Project& p) { clearPlacement(p.density); });
    });
    m_mapMode = makeButton(tr("Move the map with the mouse"), QStringLiteral("density.move_with_mouse"),
                           tr("Drag in the structure panel to move the map; Shift + drag rotates it about its "
                              "centre. Then 'Refine' to fit from there."));
    m_mapMode->setCheckable(true);
    connect(m_mapMode, &QPushButton::toggled, this, &StructurePage::mapModeChanged);
    f->addWide(m_mapMode);
    m_densityStatus = new QLabel;
    m_densityStatus->setObjectName(QStringLiteral("density.status"));
    m_densityStatus->setWordWrap(true);
    m_densityStatus->setTextInteractionFlags(Qt::TextSelectableByMouse | Qt::TextSelectableByKeyboard);
    f->addWide(m_densityStatus);
    f->addHint(tr("The map is shown on the real model: it is drawn once the window has loaded the reference "
                  "frame. The placement is stored in the project."));
}

void StructurePage::buildTags() {
    Section* s = addSection(QStringLiteral("tags"), tr("Residue number tags"));
    FormGrid* f = s->form();
    f->add(tr("Label every N-th residue"), field("visual.label_every"));
    f->add(tr("Label residues"), field("visual.label_residues", [](forms::FieldSpec& spec) {
               spec.hint = QStringLiteral("Residue numbers or ranges to tag (in addition to the interval)");
           }));
}

void StructurePage::buildAppearance() {
    Section* s = addSection(QStringLiteral("appearance"), tr("Appearance"));
    FormGrid* f = s->form();
    f->add(tr("Groups shown"), field("visual.structure_groups"));
    f->add(tr("Trace atoms (priority)"), field("visual.trace_atoms", [](forms::FieldSpec& spec) {
               spec.hint = QStringLiteral("First atom found per residue, e.g. P or C4' (CA for proteins)");
           }));
    f->addRow(tr("Options"), {check("visual.align_frames", tr("superpose frames")),
                              check("visual.show_ions", tr("show Mg²⁺")),
                              check("visual.depth_cue", tr("depth cue"))});
    f->add(tr("Colouring"), field("visual.color_mode"));
    f->addRow(tr("Base / ion colour"), {field("visual.base_color"), field("visual.ion_color")});
    f->add(tr("Chain palette"), field("visual.chain_palette"));
    f->add(tr("Proximity set"), field("visual.proximity_set"));
    f->addRow(tr("Proximity cutoff (Å) / colour"), {field("visual.proximity_cutoff"), field("visual.proximity_color")});
    f->addRow(tr("Glow"), {field("visual.glow_mode"), field("visual.glow_color")});
    f->addRow(tr("Fade frames / strength"), {field("visual.fade_frames"), field("visual.glow_strength")});
    f->addRow(tr("Line width / point size"), {field("visual.line_width"), field("visual.point_size")});
    f->addRow(tr("Theme / font size"), {field("visual.theme"), field("visual.font_size")});
    f->addRow(tr("Labels"), {check("visual.show_time_label", tr("time label")), field("visual.title")});
}

bool StructurePage::mapMode() const { return m_mapMode->isChecked(); }

void StructurePage::setCentreOnModelAvailable(bool available) {
    m_centreAvailable = available;
    m_centre->setEnabled(available);
    const QString tip = tr("Place the map's density centre on the model without rotating it (a starting point "
                           "for moving it by hand)");
    m_centre->setToolTip(available ? tip : tip + tr("\nNeeds a map file and the reference frame shown in the figure."));
}

QString StructurePage::densityStatusText() const {
    const core::DensityParams& d = context().document().project().density;
    if (d.map_path.empty()) {
        return tr("No map selected.");
    }
    if (d.transform.size() == 12) {
        return d.fit_info.empty() ? tr("Placed (not fitted).") : forms::toQString(d.fit_info);
    }
    return tr("Not placed yet: run a global search, or centre it and move it by hand.");
}

void StructurePage::documentChanged(const core::DocumentChange& change) {
    if (change.affects("density")) {
        m_densityStatus->setText(densityStatusText());
    }
    // a placement belongs to one map: another map file forgets it, in the same undo step (done
    // after the notification, which must not edit the document)
    if (m_fixingMap || change.kind != core::DocumentChange::Kind::Edit ||
        std::find(change.fields.begin(), change.fields.end(), "density.map_path") == change.fields.end() ||
        !hasPlacement(context().document().project().density)) {
        return;
    }
    const std::string label = change.label;
    const std::string mapPath = context().document().project().density.map_path;
    QPointer<StructurePage> self(this);
    QTimer::singleShot(0, this, [self, label, mapPath] {
        if (!self) {
            return;
        }
        core::ProjectDocument& doc = self->context().document();
        if (!doc.canUndo() || doc.undoLabel() != label || doc.project().density.map_path != mapPath ||
            !hasPlacement(doc.project().density)) {
            return;  // something else happened in between: leave it
        }
        self->m_fixingMap = true;
        doc.undo();
        doc.edit(label, [&mapPath](core::Project& p) {
            p.density.map_path = mapPath;
            clearPlacement(p.density);
        });
        self->m_fixingMap = false;
    });
}

}  // namespace mdbridge::app::pages
