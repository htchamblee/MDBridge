#pragma once
// The Analyses page (MDATK v0.0.4.1 gui/app.py _tab_analyses): general settings, RMSD, contacts,
// hydrogen bonds, Mg²⁺ coordination, π–π stacking, base pairs and the intermolecular energy. The
// group drop-downs (free text with suggestions) list the project's current group names and follow
// every change of the groups (MDATK _groups_changed / set_combo_entry_items).

#include <vector>

#include <QStringList>

#include "gui/pages/settings_page.hpp"

namespace mdbridge::core {
struct DocumentChange;
}

namespace mdbridge::app::forms {
class EditableComboField;
}

namespace mdbridge::app::pages {

class AnalysesPage : public SettingsPage {
    Q_OBJECT

public:
    explicit AnalysesPage(forms::FormContext& context, QWidget* parent = nullptr);

    // The group drop-downs (analyses.<contacts|hbonds|mg|stacking|energy>.group_a / group_b).
    const std::vector<forms::EditableComboField*>& groupCombos() const { return m_groupCombos; }
    // The names they offer: the project's group names, in order (blank names left out).
    QStringList groupChoices() const;

private:
    forms::EditableComboField* groupCombo(const char* analysis, const char* side);
    void updateGroupChoices();

    std::vector<forms::EditableComboField*> m_groupCombos;
    QStringList m_choices;
};

}  // namespace mdbridge::app::pages
