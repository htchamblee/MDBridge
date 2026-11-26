#pragma once
// The Residue groups page (MDATK v0.0.4.1 gui/app.py _tab_groups, with MDATK 0.0.10 / 0.1.1's group
// editing): the table of groups with colour swatches and *Pick colour…*; renaming a group renames
// it everywhere the project uses it (core::ProjectDocument::renameGroup); removing a group that is
// still used says where; *Groups from chains* adds one group per chain of the topology, and
// choosing a topology fills empty or untouched default groups with its chains (owner decision D15,
// app/chain_groups.hpp). A status line under the table says what happened (also logged).

#include <cstdint>
#include <string>
#include <vector>

#include <QString>
#include <QTimer>

#include "chain_groups.hpp"
#include "gui/pages/form_page.hpp"

class QLabel;
class QPushButton;

namespace mdbridge::core {
struct DocumentChange;
}

namespace mdbridge::app::forms {
class ListEditor;
}

namespace mdbridge::app::pages {

class GroupsPage : public FormPage {
    Q_OBJECT

public:
    explicit GroupsPage(forms::FormContext& context, QWidget* parent = nullptr);

    forms::ListEditor* groupEditor() const { return m_list; }
    QPushButton* pickColourButton() const { return m_pick; }
    QPushButton* chainsButton() const { return m_chains; }
    QLabel* statusLabel() const { return m_status; }
    // The status line: its text, its detail (tooltip; the fields a rename changed) and whether it
    // reports a problem.
    QString statusText() const;
    QString statusDetail() const;
    bool statusIsError() const { return m_statusError; }

    // Renames the group of row `row` to `text` and every reference to it (MDATK 0.1.1
    // commit_group_name): empty, invalid and taken names are refused with a message. True if renamed.
    bool commitGroupName(int row, const QString& text);
    // *Pick colour…* for the selected row (the colour dialog of forms::runColorDialog).
    void pickColour();
    // *Groups from chains* (automatic: a topology was just chosen). The chains are read on a worker
    // thread; chainsPending() until the answer arrived (then chainsDone is emitted).
    void groupsFromChains(bool automatic = false);
    bool chainsPending() const { return m_pending > 0; }

    // How long a typed topology must stay unchanged before its chains are read (MDATK: 500 ms).
    static constexpr int kTopologySettleMs = 500;

signals:
    // A line for the window's log ("WARNING: ..." for problems), as MDATK's append_log.
    void logMessage(const QString& text);
    void statusChanged(const QString& text);
    // A chain read finished (applied, refused or dropped).
    void chainsDone();

private:
    // field() with `label` (check box text, accessible name) and `hint` as tooltip when the spec has none.
    forms::FieldWidget* fieldAs(const std::string& path, const QString& label, const QString& hint = {});
    void message(const QString& text, bool error = false, const QString& detail = {});
    void onDocumentChanged(const core::DocumentChange& change);
    void topologySettled();
    void chainsReady(std::uint64_t request, const std::string& topology, const std::string& path, bool automatic,
                     const std::vector<ChainGroup>& chains, const QString& error);
    std::vector<std::string> names() const;

    forms::ListEditor* m_list = nullptr;
    QPushButton* m_pick = nullptr;
    QPushButton* m_chains = nullptr;
    QLabel* m_status = nullptr;
    bool m_statusError = false;
    std::vector<std::string> m_names;  // the group names after the last change
    bool m_renaming = false;           // a rename (or a merge of chain groups) is being stored
    QTimer m_topologyTimer;
    std::string m_offered;             // the topology whose chains were last offered
    std::uint64_t m_chainRequest = 0;  // the newest read; older answers are dropped
    std::uint64_t m_chainsReading = 0; // the read the status line says is running
    QString m_readingText;             // "Reading the chains of <file> …"
    QString m_readingName;
    int m_pending = 0;
};

}  // namespace mdbridge::app::pages
