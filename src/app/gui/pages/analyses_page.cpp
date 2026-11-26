#include "gui/pages/analyses_page.hpp"

#include <string>

#include <QLabel>
#include <QTableView>

#include "gui/forms/field_widgets.hpp"
#include "gui/forms/list_editor.hpp"
#include "mdbridge/core/project_document.hpp"

namespace mdbridge::app::pages {

AnalysesPage::AnalysesPage(forms::FormContext& context, QWidget* parent)
    : SettingsPage(context, QStringLiteral("analyses"), QStringLiteral("Analyses"), parent) {
    // ---- General
    Section* gen = addSection(QStringLiteral("general"), tr("General"));
    gen->addWide(field("analyses.use_pbc",
                       tr("Use periodic (minimum-image) distances (not needed after cpptraj imaging)")));
    gen->add(tr("Proximity search radius (Å)"),
             field("analyses.proximity_search_radius", tr("Proximity search radius (Å)"),
                   tr("Largest proximity cutoff usable at render time without recomputing")),
             tr("Largest proximity cutoff usable at render time without recomputing"));

    // ---- RMSD
    Section* r = addSection(QStringLiteral("rmsd"), tr("RMSD"));
    r->addWide(field("analyses.rmsd.enabled", tr("Compute RMSD")));
    forms::FieldWidget* targets = field("analyses.rmsd.targets", tr("Targets"));
    if (auto* list = qobject_cast<forms::ListEditor*>(targets)) {
        list->view()->setMinimumHeight(90);
    }
    r->add(tr("Targets"), targets);
    const QString refHint = tr("Index into the analysed frames");
    r->addRow(tr("Reference frame"),
              {field("analyses.rmsd.reference_frame", tr("Reference frame"), refHint),
               field("analyses.rmsd.superposition", tr("superpose")),
               field("analyses.rmsd.mass_weighted", tr("mass weighted"))},
              refHint);

    // ---- Pairwise contacts
    Section* c = addSection(QStringLiteral("contacts"), tr("Pairwise contacts"));
    c->addWide(field("analyses.contacts.enabled", tr("Compute heavy-atom contacts")));
    c->addRow(tr("Groups"), {groupCombo("contacts", "group_a"), rowLabel(tr("vs")), groupCombo("contacts", "group_b")});
    c->addRow(tr("Cutoff (Å)"), {field("analyses.contacts.cutoff", tr("Cutoff (Å)")),
                                 field("analyses.contacts.heavy_only", tr("heavy atoms only"))});

    // ---- Hydrogen bonds
    Section* h = addSection(QStringLiteral("hbonds"), tr("Hydrogen bonds"));
    h->addWide(field("analyses.hbonds.enabled", tr("Count inter-group hydrogen bonds")));
    h->addRow(tr("Groups"), {groupCombo("hbonds", "group_a"), rowLabel(tr("vs")), groupCombo("hbonds", "group_b")});
    h->addRow(tr("D–A distance (Å) / angle"),
              {field("analyses.hbonds.distance_cutoff", tr("D–A distance (Å)")),
               field("analyses.hbonds.use_angle", tr("use D–H···A angle ≥")),
               field("analyses.hbonds.angle_cutoff", tr("D–H···A angle"))});
    h->addRow(tr("Donor / acceptor elements"), {field("analyses.hbonds.donor_elements", tr("Donor elements")),
                                                field("analyses.hbonds.acceptor_elements", tr("Acceptor elements"))});

    // ---- O–Mg–O contacts and Mg occupancy
    Section* m = addSection(QStringLiteral("mg"), tr("O–Mg–O contacts and Mg occupancy"));
    m->addWide(field("analyses.mg.enabled", tr("Analyse Mg²⁺ coordination")));
    m->add(tr("Mg selection"), field("analyses.mg.mg_selection", tr("Mg selection")));
    const QString ligandHint = tr("Oxygens of these residue ranges count for Mg–O and O–Mg–O");
    m->addRow(tr("Ligand selection / elements"),
              {field("analyses.mg.ligand_selection", tr("Ligand selection"), ligandHint),
               field("analyses.mg.ligand_elements", tr("Ligand elements"), ligandHint)},
              ligandHint);
    m->addRow(tr("Bridging groups"), {groupCombo("mg", "group_a"), rowLabel(tr("vs")), groupCombo("mg", "group_b")});
    m->addRow(tr("Inner / outer sphere (Å)"), {field("analyses.mg.inner_cutoff", tr("Inner sphere (Å)")),
                                               field("analyses.mg.outer_cutoff", tr("Outer sphere (Å)"))});
    m->add(tr("Bridges counted from"), field("analyses.mg.bridge_shell", tr("Bridges counted from")));
    m->addRow(tr("Coordination geometry"),
              {field("analyses.mg.geometry", tr("require octahedral O–Mg–O angles ±")),
               field("analyses.mg.angle_tolerance", tr("Angle tolerance"))});
    m->addRow(tr("Occupancy residues / atoms"),
              {field("analyses.mg.occupancy_selection", tr("Occupancy residues")),
               field("analyses.mg.occupancy_elements", tr("Occupancy atoms"))});

    // ---- π–π stacking
    Section* s = addSection(QStringLiteral("stacking"), tr("π–π stacking"));
    s->addWide(field("analyses.stacking.enabled", tr("Detect π–π stacking")));
    s->add(tr("Selection"), field("analyses.stacking.selection", tr("Selection")));
    s->add(tr("Scope"), field("analyses.stacking.scope", tr("Scope")));
    s->addRow(tr("Groups (A vs B)"),
              {groupCombo("stacking", "group_a"), rowLabel(tr("vs")), groupCombo("stacking", "group_b")},
              tr("Used when the scope is 'group A vs group B'"));
    const QString centroidHint = tr("Ring-centroid distance and angle between ring normals");
    s->addRow(tr("Centroid distance (Å) / angle (°)"),
              {field("analyses.stacking.distance_cutoff", tr("Centroid distance (Å)"), centroidHint),
               field("analyses.stacking.angle_cutoff", tr("Angle between ring normals (°)"), centroidHint)},
              centroidHint);
    const QString offsetHint = tr("Lateral centroid offset and RMS deviation of the ring atoms from their plane");
    s->addRow(tr("Offset / planarity (Å)"),
              {field("analyses.stacking.offset_cutoff", tr("Offset (Å)"), offsetHint),
               field("analyses.stacking.planarity_cutoff", tr("Planarity (Å)"), offsetHint)},
              offsetHint);

    // ---- Base pairs (secondary structure)
    Section* bp = addSection(QStringLiteral("basepairs"), tr("Base pairs (secondary structure)"));
    bp->addWide(hintLabel(tr("Used by the secondary-structure panel (and PDB annotation): pairs of the "
                             "structure-view residues, re-detected every n frames.")));
    bp->add(tr("Engine"), field("analyses.basepairs.engine", tr("Engine")));
    const QString strideHint = tr("Pairs are re-detected every n-th analysed frame; the cutoff is the built-in "
                                  "engine's largest donor-acceptor distance");
    bp->addRow(tr("Every n frames / cutoff (Å)"),
               {field("analyses.basepairs.stride", tr("Every n frames"), strideHint),
                field("analyses.basepairs.cutoff", tr("Cutoff (Å)"), strideHint)},
               strideHint);
    bp->addWide(field("analyses.basepairs.wobble", tr("Built-in: count G-U wobble pairs")));
    bp->addWide(field("analyses.basepairs.noncanonical", tr("Barnaba: keep non-canonical pairs")));

    // ---- Intermolecular energy
    Section* e = addSection(QStringLiteral("energy"), tr("Intermolecular energy"));
    e->addWide(field("analyses.energy.enabled", tr("Compute intermolecular interaction energy and net force")));
    e->addRow(tr("Groups"), {groupCombo("energy", "group_a"), rowLabel(tr("vs")), groupCombo("energy", "group_b")});
    e->addRow(tr("Cutoff (Å) / dielectric"), {field("analyses.energy.cutoff", tr("Cutoff (Å)")),
                                              field("analyses.energy.dielectric", tr("Dielectric"))});
    e->addWide(field("analyses.energy.distance_dependent", tr("Distance-dependent dielectric (εr)")));
    e->addWide(field("analyses.energy.use_parm7", tr("Amber LJ parameters from the parm7 when available")));
    e->addRow(tr("Exclude pairs within"),
              {field("analyses.energy.exclude_bonded", tr("Exclude pairs within (bonds)")), rowLabel(tr("bonds"))});
    e->add(tr("Target value (kcal/mol)"), field("analyses.energy.target_value", tr("Target value (kcal/mol)")),
           tr("Drawn as a reference line on charts and panels"));
    finish();

    updateGroupChoices();
    connect(&context, &forms::FormContext::documentChanged, this, [this](const core::DocumentChange& change) {
        if (change.affects("groups")) {
            updateGroupChoices();
        }
    });
}

forms::EditableComboField* AnalysesPage::groupCombo(const char* analysis, const char* side) {
    const std::string path = std::string("analyses.") + analysis + "." + side;
    auto* combo = qobject_cast<forms::EditableComboField*>(
        field(path, std::string(side) == "group_a" ? tr("Group A") : tr("Group B")));
    m_groupCombos.push_back(combo);
    return combo;
}

QStringList AnalysesPage::groupChoices() const {
    QStringList out;
    for (const core::ResidueGroup& g : context().document().project().groups) {
        const QString name = forms::toQString(g.name).trimmed();
        if (!name.isEmpty()) out.append(name);
    }
    return out;
}

void AnalysesPage::updateGroupChoices() {
    const QStringList names = groupChoices();
    if (names == m_choices) {
        return;
    }
    m_choices = names;
    for (forms::EditableComboField* combo : m_groupCombos) {
        combo->setChoices(names);
    }
}

}  // namespace mdbridge::app::pages
