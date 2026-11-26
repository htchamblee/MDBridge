#include "gui/pages/analyses_page.hpp"

#include <string>

#include <QLabel>
#include <QTableView>

#include "gui/forms/field_widgets.hpp"
#include "gui/forms/list_editor.hpp"
#include "mdbridge/core/project_document.hpp"

namespace mdbridge::app::pages {

AnalysesPage::AnalysesPage(forms::FormContext& context, QWidget* parent)
    : FormPage(context, QStringLiteral("analyses"), QStringLiteral("Analyses"), parent) {
    // ---- General
    Section* gen = addSection(QStringLiteral("general"), tr("General"));
    gen->form()->addWide(fieldAs("analyses.use_pbc",
                       tr("Use periodic (minimum-image) distances (not needed after cpptraj imaging)")));
    gen->form()->add(tr("Proximity search radius (Å)"),
             fieldAs("analyses.proximity_search_radius", tr("Proximity search radius (Å)"),
                   tr("Largest proximity cutoff usable at render time without recomputing")),
             tr("Largest proximity cutoff usable at render time without recomputing"));

    // ---- RMSD
    Section* r = addSection(QStringLiteral("rmsd"), tr("RMSD"));
    r->form()->addWide(fieldAs("analyses.rmsd.enabled", tr("Compute RMSD")));
    forms::FieldWidget* targets = fieldAs("analyses.rmsd.targets", tr("Targets"));
    if (auto* list = qobject_cast<forms::ListEditor*>(targets)) {
        list->view()->setMinimumHeight(90);
    }
    r->form()->add(tr("Targets"), targets);
    const QString refHint = tr("Index into the analysed frames");
    r->form()->addRow(tr("Reference frame"),
              {fieldAs("analyses.rmsd.reference_frame", tr("Reference frame"), refHint),
               fieldAs("analyses.rmsd.superposition", tr("superpose")),
               fieldAs("analyses.rmsd.mass_weighted", tr("mass weighted"))},
              refHint);

    // ---- Pairwise contacts
    Section* c = addSection(QStringLiteral("contacts"), tr("Pairwise contacts"));
    c->form()->addWide(fieldAs("analyses.contacts.enabled", tr("Compute heavy-atom contacts")));
    c->form()->addRow(tr("Groups"), {groupCombo("contacts", "group_a"), new QLabel(tr("vs")), groupCombo("contacts", "group_b")});
    c->form()->addRow(tr("Cutoff (Å)"), {fieldAs("analyses.contacts.cutoff", tr("Cutoff (Å)")),
                                 fieldAs("analyses.contacts.heavy_only", tr("heavy atoms only"))});

    // ---- Hydrogen bonds
    Section* h = addSection(QStringLiteral("hbonds"), tr("Hydrogen bonds"));
    h->form()->addWide(fieldAs("analyses.hbonds.enabled", tr("Count inter-group hydrogen bonds")));
    h->form()->addRow(tr("Groups"), {groupCombo("hbonds", "group_a"), new QLabel(tr("vs")), groupCombo("hbonds", "group_b")});
    h->form()->addRow(tr("D–A distance (Å) / angle"),
              {fieldAs("analyses.hbonds.distance_cutoff", tr("D–A distance (Å)")),
               fieldAs("analyses.hbonds.use_angle", tr("use D–H···A angle ≥")),
               fieldAs("analyses.hbonds.angle_cutoff", tr("D–H···A angle"))});
    h->form()->addRow(tr("Donor / acceptor elements"), {fieldAs("analyses.hbonds.donor_elements", tr("Donor elements")),
                                                fieldAs("analyses.hbonds.acceptor_elements", tr("Acceptor elements"))});

    // ---- O–Mg–O contacts and Mg occupancy
    Section* m = addSection(QStringLiteral("mg"), tr("O–Mg–O contacts and Mg occupancy"));
    m->form()->addWide(fieldAs("analyses.mg.enabled", tr("Analyse Mg²⁺ coordination")));
    m->form()->add(tr("Mg selection"), fieldAs("analyses.mg.mg_selection", tr("Mg selection")));
    const QString ligandHint = tr("Oxygens of these residue ranges count for Mg–O and O–Mg–O");
    m->form()->addRow(tr("Ligand selection / elements"),
              {fieldAs("analyses.mg.ligand_selection", tr("Ligand selection"), ligandHint),
               fieldAs("analyses.mg.ligand_elements", tr("Ligand elements"), ligandHint)},
              ligandHint);
    m->form()->addRow(tr("Bridging groups"), {groupCombo("mg", "group_a"), new QLabel(tr("vs")), groupCombo("mg", "group_b")});
    m->form()->addRow(tr("Inner / outer sphere (Å)"), {fieldAs("analyses.mg.inner_cutoff", tr("Inner sphere (Å)")),
                                               fieldAs("analyses.mg.outer_cutoff", tr("Outer sphere (Å)"))});
    m->form()->add(tr("Bridges counted from"), fieldAs("analyses.mg.bridge_shell", tr("Bridges counted from")));
    m->form()->addRow(tr("Coordination geometry"),
              {fieldAs("analyses.mg.geometry", tr("require octahedral O–Mg–O angles ±")),
               fieldAs("analyses.mg.angle_tolerance", tr("Angle tolerance"))});
    m->form()->addRow(tr("Occupancy residues / atoms"),
              {fieldAs("analyses.mg.occupancy_selection", tr("Occupancy residues")),
               fieldAs("analyses.mg.occupancy_elements", tr("Occupancy atoms"))});

    // ---- π–π stacking
    Section* s = addSection(QStringLiteral("stacking"), tr("π–π stacking"));
    s->form()->addWide(fieldAs("analyses.stacking.enabled", tr("Detect π–π stacking")));
    s->form()->add(tr("Selection"), fieldAs("analyses.stacking.selection", tr("Selection")));
    s->form()->add(tr("Scope"), fieldAs("analyses.stacking.scope", tr("Scope")));
    s->form()->addRow(tr("Groups (A vs B)"),
              {groupCombo("stacking", "group_a"), new QLabel(tr("vs")), groupCombo("stacking", "group_b")},
              tr("Used when the scope is 'group A vs group B'"));
    const QString centroidHint = tr("Ring-centroid distance and angle between ring normals");
    s->form()->addRow(tr("Centroid distance (Å) / angle (°)"),
              {fieldAs("analyses.stacking.distance_cutoff", tr("Centroid distance (Å)"), centroidHint),
               fieldAs("analyses.stacking.angle_cutoff", tr("Angle between ring normals (°)"), centroidHint)},
              centroidHint);
    const QString offsetHint = tr("Lateral centroid offset and RMS deviation of the ring atoms from their plane");
    s->form()->addRow(tr("Offset / planarity (Å)"),
              {fieldAs("analyses.stacking.offset_cutoff", tr("Offset (Å)"), offsetHint),
               fieldAs("analyses.stacking.planarity_cutoff", tr("Planarity (Å)"), offsetHint)},
              offsetHint);

    // ---- Base pairs (secondary structure)
    Section* bp = addSection(QStringLiteral("basepairs"), tr("Base pairs (secondary structure)"));
    bp->form()->addWide(makeHint(tr("Used by the secondary-structure panel (and PDB annotation): pairs of the "
                             "structure-view residues, re-detected every n frames.")));
    bp->form()->add(tr("Engine"), fieldAs("analyses.basepairs.engine", tr("Engine")));
    const QString strideHint = tr("Pairs are re-detected every n-th analysed frame; the cutoff is the built-in "
                                  "engine's largest donor-acceptor distance");
    bp->form()->addRow(tr("Every n frames / cutoff (Å)"),
               {fieldAs("analyses.basepairs.stride", tr("Every n frames"), strideHint),
                fieldAs("analyses.basepairs.cutoff", tr("Cutoff (Å)"), strideHint)},
               strideHint);
    bp->form()->addWide(fieldAs("analyses.basepairs.wobble", tr("Built-in: count G-U wobble pairs")));
    bp->form()->addWide(fieldAs("analyses.basepairs.noncanonical", tr("Barnaba: keep non-canonical pairs")));

    // ---- Intermolecular energy
    Section* e = addSection(QStringLiteral("energy"), tr("Intermolecular energy"));
    e->form()->addWide(fieldAs("analyses.energy.enabled", tr("Compute intermolecular interaction energy and net force")));
    e->form()->addRow(tr("Groups"), {groupCombo("energy", "group_a"), new QLabel(tr("vs")), groupCombo("energy", "group_b")});
    e->form()->addRow(tr("Cutoff (Å) / dielectric"), {fieldAs("analyses.energy.cutoff", tr("Cutoff (Å)")),
                                              fieldAs("analyses.energy.dielectric", tr("Dielectric"))});
    e->form()->addWide(fieldAs("analyses.energy.distance_dependent", tr("Distance-dependent dielectric (εr)")));
    e->form()->addWide(fieldAs("analyses.energy.use_parm7", tr("Amber LJ parameters from the parm7 when available")));
    e->form()->addRow(tr("Exclude pairs within"),
              {fieldAs("analyses.energy.exclude_bonded", tr("Exclude pairs within (bonds)")), new QLabel(tr("bonds"))});
    e->form()->add(tr("Target value (kcal/mol)"), fieldAs("analyses.energy.target_value", tr("Target value (kcal/mol)")),
           tr("Drawn as a reference line on charts and panels"));

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
        fieldAs(path, std::string(side) == "group_a" ? tr("Group A") : tr("Group B")));
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

forms::FieldWidget* AnalysesPage::fieldAs(const std::string& path, const QString& label, const QString& hint) {
    return field(path, [&](forms::FieldSpec& spec) {
        spec.label = label;
        if (spec.hint.isEmpty()) {
            spec.hint = hint;
        }
    });
}

}  // namespace mdbridge::app::pages
