#pragma once
// The "Layout & panels" settings page (milestone 5 W2-A; MDATK v0.0.4.1 gui/app.py:1109-1334
// _tab_layout, _slot_editor, _read_slots): the grid preset and one editor per slot of it (panel
// type, data key with its label, the residue groups of the panel's analysis with "(analysis
// default)", more keys, the title, the options as JSON), and the real-data still of a frame.
//
// Bound to the document: an editor shows the panel of its slot and stores each change as one undo
// step. The editors work as MDATK's: the data drop-down lists the keys of the panel type's kind
// (expected_keys) labelled "<label>   [<key>]"; the groups make the data key an analysis instance
// ("contacts@Aptamer+Mg/atom_pairs", the analysis's own groups filling a slot left at its
// default); the title follows the data until a title of one's own is typed (a new type always gets
// its own title; the default titles count as automatic); a new sequence or secondary-structure
// panel gets MDATK's default options. Deliberate differences (docs/DIFFERENCES.md):
//
// - AP-07: a slot without a panel stays empty in the file; it becomes a panel when its editor is
//   changed or a panel is dragged onto it (and so do the empty slots before it, which hold what
//   their editors show).
// - AP-08: a preset with fewer slots keeps the panels beyond them; a larger preset shows them again.
// - AP-09 (MDATK 0.1.0): a swap or a preset change while a slot's options are no JSON object is
//   refused ("Panels not swapped; fix this first: ...", "Grid preset not changed; fix this first:
//   ..." with the preset put back).
// - AP-10: options that do not parse are a located form error ("Layout & panels › Panel 1
//   (Structure) › Options: '{bad' is not valid JSON: ..."); the document keeps the stored options.
// - Automatic titles follow a group rename (MDATK 0.1.1 _adopt_project): retitleAutomaticPanels.

#include <functional>
#include <memory>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include <QString>

#include "gui/pages/form_page.hpp"
#include "mdbridge/core/project.hpp"

class QComboBox;
class QLineEdit;
class QPushButton;
class QSpinBox;
class QVBoxLayout;

namespace mdbridge::app::pages {

class SlotTextField;

// MDATK 0.1.1's automatic titles after a rename (_adopt_project): a panel whose title is its data's
// label in `before` (or the automatic or a default title of a panel whose keys are among `changed`)
// gets its data's label in `project`. Returns the paths of the titles changed
// ("layout.panels.2.title").
std::vector<std::string> retitleAutomaticPanels(core::Project& project, const core::Project& before,
                                                const std::vector<std::string>& changed);

class LayoutPage : public FormPage {
    Q_OBJECT

public:
    explicit LayoutPage(forms::FormContext& context, QWidget* parent = nullptr);
    ~LayoutPage() override;

    // What slot editor `k` shows (the GUI tests' snapshots: tests/data/README.md "gui").
    struct SlotState {
        QString frame;                // "Slot 1 (large)"
        QString type;
        std::optional<QString> data;  // nullopt: the drop-down has no item
        std::vector<std::pair<QString, QString>> choices;  // id, label
        QString groupA, groupB;       // "" = (analysis default)
        bool groupAEnabled = false, groupBEnabled = false;
        QString more, title, options;
    };
    int slotCount() const { return static_cast<int>(m_slots.size()); }
    SlotState slotState(int k) const;

    // A panel dragged onto another slot (the canvas): swapped as one undo step and logged
    // ("Swapped panels 2 and 3"); refused while a slot's options do not parse (AP-09). A panel
    // swapped into an empty slot makes it (and the slots before it) hold panels (AP-07).
    bool swapPanels(int from, int to);
    // The located message of the first slot whose options do not parse ("" when all do).
    QString optionsProblem() const;

    QComboBox* presetBox() const { return m_preset; }
    QSpinBox* stillFrame() const { return m_stillFrame; }
    QPushButton* renderStillButton() const { return m_renderStill; }

signals:
    // A line for the window's log (MDATK's append_log: refused swaps and presets, swaps).
    void logLine(const QString& text);
    // The preview could not be updated: `problem` (logged "Preview not updated: <problem>").
    void previewNotUpdated(const QString& problem);
    // Redraw the preview now (a slot's text was committed, MDATK's activate -> refresh_canvas).
    void previewRequested();
    // "Render still" of frame `frame` (wired by the window's actions).
    void renderStillRequested(int frame);

private:
    struct Slot;
    friend class SlotTextField;

    void rebuildSlots();
    void fixTabOrder();
    void refreshSlots(bool force);
    void loadSlot(Slot& slot, const core::PanelSpec& spec, bool keepText);
    void fillData(Slot& slot, const std::string& type, const std::string& select);
    void updateGroupWidgets(Slot& slot);
    void fillGroups(Slot& slot);
    std::string primaryKey(const Slot& slot) const;
    std::string autoTitle(const Slot& slot) const;
    core::PanelSpec defaultSpec() const;
    core::PanelSpec editorSpec(const Slot& slot) const;  // what the editor shows, as stored
    bool editorMatches(const Slot& slot, const core::PanelSpec& spec) const;
    void onPreset(int index);
    void onType(Slot& slot);
    void onData(Slot& slot, bool newType);
    void commitText(Slot& slot, const QString& field);
    // Stores `change` of slot k's panel as one step (filling the empty slots up to k, AP-07).
    void storeSlot(int k, const QString& label, const std::function<void(core::PanelSpec&)>& change);
    void documentChanged(const core::DocumentChange& change);

    QComboBox* m_preset = nullptr;
    QVBoxLayout* m_slotBox = nullptr;
    QSpinBox* m_stillFrame = nullptr;
    QPushButton* m_renderStill = nullptr;
    std::vector<std::unique_ptr<Slot>> m_slots;
    int m_loading = 0;
};

}  // namespace mdbridge::app::pages
