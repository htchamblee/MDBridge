#include "gui/pages/groups_page.hpp"
#include "gui/pages/layout_page.hpp"

#include <algorithm>
#include <exception>
#include <filesystem>
#include <memory>
#include <stdexcept>
#include <utility>

#include <QApplication>
#include <QColor>
#include <QHBoxLayout>
#include <QLabel>
#include <QPalette>
#include <QPointer>
#include <QPushButton>
#include <QStyle>
#include <QTableView>
#include <QThreadPool>

#include "gui/forms/field_widgets.hpp"
#include "gui/forms/list_editor.hpp"
#include "mdbridge/core/errors.hpp"
#include "mdbridge/core/paths.hpp"
#include "mdbridge/core/project_document.hpp"

namespace mdbridge::app::pages {

namespace {

std::string strip(const std::string& text) {
    const std::size_t a = text.find_first_not_of(" \t\r\n");
    if (a == std::string::npos) return {};
    return text.substr(a, text.find_last_not_of(" \t\r\n") - a + 1);
}

QString joined(const std::vector<std::string>& items, const char* separator = ", ") {
    QStringList out;
    for (const std::string& s : items) out.append(forms::toQString(s));
    return out.join(QLatin1String(separator));
}

QString fileName(const std::string& path) {
    return forms::toQString(core::pathToUtf8(core::pathFromUtf8(path).filename()));
}

QString readingText(const std::string& path) {
    return QObject::tr("Reading the chains of %1 …").arg(fileName(path));
}

}  // namespace

GroupsPage::GroupsPage(forms::FormContext& context, QWidget* parent)
    : FormPage(context, QStringLiteral("groups"), QStringLiteral("Residue groups"), parent) {
    Section* s = addSection(QStringLiteral("groups"), tr("Residue groups"));
    m_list = qobject_cast<forms::ListEditor*>(fieldAs("groups", tr("Residue groups")));
    m_list->view()->setMinimumHeight(180);
    m_pick = makeButton(tr("Pick colour…"), QStringLiteral("groups:pick_colour"), tr("Pick the selected group's colour"));
    m_pick->setIcon(style()->standardIcon(QStyle::SP_DialogApplyButton));
    m_list->addButton(m_pick);
    m_list->setCellEditHandler([this](int row, const std::string& key, const QString& text) {
        if (key != "name") {
            return false;
        }
        commitGroupName(row, text);
        return true;  // stored, or refused with a message: never stored as plain text
    });
    s->form()->addWide(m_list);

    auto* chainsRow = new QWidget;
    auto* row = new QHBoxLayout(chainsRow);
    row->setContentsMargins(0, 0, 0, 0);
    m_chains = makeButton(tr("Groups from chains"), QStringLiteral("groups:from_chains"),
                          tr("Add one group per chain of the topology (water and ions left out); existing groups "
                             "are kept"));
    m_chains->setIcon(style()->standardIcon(QStyle::SP_FileDialogNewFolder));
    row->addWidget(m_chains);
    row->addStretch(1);
    s->form()->addWide(chainsRow);
    m_status = new QLabel;
    m_status->setObjectName(QStringLiteral("groups:status"));
    m_status->setWordWrap(true);
    m_status->setTextInteractionFlags(Qt::TextSelectableByMouse | Qt::TextSelectableByKeyboard);
    s->form()->addWide(m_status);
    s->form()->addWide(makeHint(tr(
        "Residues accept ranges such as “1-40, 45, 50-60” or chain/segment-qualified “A:1-40”. A "
        "selection is a raw MDAnalysis selection string. Anywhere a selection is expected you can refer to a group "
        "as group:<Name> (e.g. “group:Aptamer and name P”). Group names must not contain spaces. Renaming "
        "a group (Enter or leaving the cell) renames it everywhere the project uses it. Choosing a topology fills "
        "empty or untouched default groups with one group per chain.")));

    m_names = names();
    m_topologyTimer.setSingleShot(true);
    m_topologyTimer.setInterval(kTopologySettleMs);
    connect(&m_topologyTimer, &QTimer::timeout, this, &GroupsPage::topologySettled);
    connect(m_pick, &QPushButton::clicked, this, &GroupsPage::pickColour);
    connect(m_chains, &QPushButton::clicked, this, [this] { groupsFromChains(false); });
    connect(&context, &forms::FormContext::documentChanged, this, &GroupsPage::onDocumentChanged);
}

QString GroupsPage::statusText() const { return m_status->text(); }
QString GroupsPage::statusDetail() const { return m_status->toolTip(); }

std::vector<std::string> GroupsPage::names() const {
    std::vector<std::string> out;
    for (const core::ResidueGroup& g : context().document().project().groups) {
        if (!strip(g.name).empty()) out.push_back(strip(g.name));
    }
    return out;
}

void GroupsPage::message(const QString& text, bool error, const QString& detail) {
    m_status->setText(text);
    m_status->setToolTip(detail);
    m_statusError = error;
    QPalette palette = m_status->palette();
    palette.setColor(QPalette::WindowText, error ? QColor(0xc0, 0x1c, 0x28) : QApplication::palette().color(QPalette::WindowText));
    m_status->setPalette(palette);
    if (!text.isEmpty()) {
        emit logMessage((error ? QStringLiteral("WARNING: ") : QString()) + text +
                        (detail.isEmpty() ? QString() : QStringLiteral("\n  ") + detail));
    }
    emit statusChanged(text);
}

bool GroupsPage::commitGroupName(int row, const QString& text) {
    core::ProjectDocument& doc = context().document();
    const std::vector<core::ResidueGroup>& groups = doc.project().groups;
    if (row < 0 || static_cast<std::size_t>(row) >= groups.size()) {
        return false;
    }
    const std::string old = strip(groups[static_cast<std::size_t>(row)].name);
    const std::string name = strip(forms::toStdString(text));
    if (name == old) {
        return false;
    }
    std::vector<std::string> others;
    for (std::size_t k = 0; k < groups.size(); ++k) {
        if (k != static_cast<std::size_t>(row)) others.push_back(strip(groups[k].name));
    }
    if (const std::string problem = core::groupNameProblem(name, others); !problem.empty()) {
        message(tr("Group %1 not renamed: %2").arg(old.empty() ? tr("(unnamed)") : forms::toQString(old),
                                                   forms::toQString(problem)),
                true);
        return false;
    }
    m_renaming = true;
    const auto done = [this] { m_renaming = false; };
    const int same = static_cast<int>(std::count(others.begin(), others.end(), old));
    try {
        if (old.empty() || same > 0) {
            // not a group yet, or one of several of that name: no references to follow
            doc.set("groups." + std::to_string(row) + ".name", core::Value(name));
            done();
            message(old.empty() ? tr("Named the new group %1").arg(forms::toQString(name))
                                : tr("Renamed %1 → %2").arg(forms::toQString(old), forms::toQString(name)));
            return true;
        }
        // one step: the rename and the automatic panel titles that follow it (MDATK 0.1.1)
        const core::Project before = doc.project();
        std::vector<std::string> refs;
        doc.beginGroup("Rename group " + old);
        try {
            refs = doc.renameGroup(old, name);
            core::Project next = doc.project();
            const std::vector<std::string> titles = retitleAutomaticPanels(next, before, refs);
            if (!titles.empty()) {
                doc.edit("Retitle panels", [&](core::Project& p) { p = next; });
                refs.insert(refs.end(), titles.begin(), titles.end());
            }
        } catch (...) {
            doc.cancelGroup();
            throw;
        }
        doc.endGroup();
        done();
        message((refs.size() == 1 ? tr("Renamed %1 → %2 in %3 place") : tr("Renamed %1 → %2 in %3 places"))
                    .arg(forms::toQString(old), forms::toQString(name), QString::number(refs.size())),
                false, joined(refs));
        return true;
    } catch (const std::exception& e) {
        done();
        message(tr("Group %1 not renamed; fix these settings first: %2").arg(forms::toQString(old), QString::fromUtf8(e.what())),
                true);
        return false;
    }
}

void GroupsPage::pickColour() {
    const int row = m_list->selectedRow();
    core::ProjectDocument& doc = context().document();
    if (row < 0 || static_cast<std::size_t>(row) >= doc.project().groups.size()) {
        message(tr("Select a group first."), true);
        return;
    }
    const std::string path = "groups." + std::to_string(row) + ".color";
    const QColor chosen = forms::runColorDialog(window(), forms::matplotlibColor(forms::toQString(doc.text(path))),
                                                tr("Group colour"));
    if (chosen.isValid()) {
        doc.set(path, core::Value(forms::toStdString(chosen.name(QColor::HexRgb))), "Pick group colour");
    }
}

void GroupsPage::onDocumentChanged(const core::DocumentChange& change) {
    using Kind = core::DocumentChange::Kind;
    if (change.kind == Kind::Replaced) {
        // a project loaded (New, Open): pending reads and offers belong to the one before
        ++m_chainRequest;
        m_topologyTimer.stop();
        m_offered.clear();
        m_names = names();
        message({});
        return;
    }
    if (change.affects("groups")) {
        const std::vector<std::string> now = names();
        std::vector<std::string> gone;
        for (const std::string& n : m_names) {
            if (std::find(now.begin(), now.end(), n) == now.end()) gone.push_back(n);
        }
        m_names = now;
        if (!gone.empty() && change.kind == Kind::Edit && !m_renaming) {
            // say where removed groups are still used: those settings now name a group that is gone
            QStringList uses;
            for (const std::string& n : gone) {
                const std::vector<std::string> refs = context().document().groupReferences(n);
                if (!refs.empty()) uses.append(tr("%1 in %2").arg(forms::toQString(n), joined(refs)));
            }
            if (!uses.isEmpty()) {
                message(tr("Removed group still in use, change these settings: ") + uses.join(QStringLiteral("; ")), true);
            }
        }
    }
    if (change.kind == Kind::Edit && change.affects("inputs.topology")) {
        m_topologyTimer.start();  // offer the chains once the name stops changing
    }
}

void GroupsPage::topologySettled() {
    const std::string text = strip(context().document().text("inputs.topology"));
    if (text.empty() || text == m_offered) {
        return;
    }
    std::string path;
    try {
        path = context().document().resolvedProject().inputs.topology;
    } catch (const std::exception&) {
        return;
    }
    std::error_code ec;
    if (std::filesystem::is_regular_file(core::pathFromUtf8(path), ec)) {
        m_offered = text;
        groupsFromChains(true);
    }
}

void GroupsPage::groupsFromChains(bool automatic) {
    core::ProjectDocument& doc = context().document();
    std::string path;
    try {
        path = doc.resolvedProject().inputs.topology;
    } catch (const std::exception&) {
        path.clear();
    }
    std::error_code ec;
    if (strip(doc.text("inputs.topology")).empty() || !std::filesystem::is_regular_file(core::pathFromUtf8(path), ec)) {
        if (!automatic) {
            message(tr("Choose a topology file on the Input page first."), true);
        }
        return;
    }
    const std::uint64_t request = ++m_chainRequest;
    const std::vector<core::ResidueGroup> existing = doc.project().groups;
    const std::string topology = doc.text("inputs.topology");
    if (!automatic) {
        m_chainsReading = request;
        m_readingText = readingText(path);
        m_readingName = fileName(path);
        m_status->setText(m_readingText);
        m_status->setToolTip({});
        emit statusChanged(m_status->text());
    }
    ++m_pending;
    QPointer<GroupsPage> self(this);
    QThreadPool::globalInstance()->start([self, request, topology, automatic, path, existing] {
        std::vector<ChainGroup> chains;
        QString error;
        try {
            chains = topologyChainGroups(core::pathFromUtf8(path), existing);
        } catch (const std::exception& e) {  // every reader problem is reported
            error = QString::fromUtf8(e.what());
        }
        // back on the GUI thread; dropped when the page is gone
        QMetaObject::invokeMethod(
            qApp,
            [self, request, topology, path, automatic, chains = std::move(chains), error] {
                if (self) {
                    self->chainsReady(request, topology, path, automatic, chains, error);
                }
            },
            Qt::QueuedConnection);
    });
}

void GroupsPage::chainsReady(std::uint64_t request, const std::string& topology, const std::string& path, bool automatic,
                             const std::vector<ChainGroup>& chains, const QString& error) {
    --m_pending;
    core::ProjectDocument& doc = context().document();
    const QString name = fileName(path);
    // superseded by a newer request or by loading a project, or the topology changed meanwhile
    if (request != m_chainRequest || topology != doc.text("inputs.topology")) {
        if (request == m_chainsReading && m_status->text() == m_readingText) {
            message(tr("Groups not changed: the topology changed while the chains of %1 were read").arg(m_readingName));
        }
        emit chainsDone();
        return;
    }
    if (!error.isEmpty()) {
        message(tr("Could not read the chains of %1: %2").arg(name, error), true);
        emit chainsDone();
        return;
    }
    core::Project merged = doc.project();
    const ChainMerge merge = mergeChainGroups(merged, chains, automatic);
    if (merge.changed()) {
        m_renaming = true;  // replaced defaults are renamed, not removed
        try {
            doc.edit("Groups from chains", [&](core::Project& p) { p = merged; });
        } catch (...) {
            m_renaming = false;
            throw;
        }
        m_renaming = false;
    }
    if (merge.changed() || !merge.pending.empty() || !automatic) {
        message(forms::toQString(merge.summary(forms::toStdString(name))), false, joined(merge.references));
    }
    emit chainsDone();
}

forms::FieldWidget* GroupsPage::fieldAs(const std::string& path, const QString& label, const QString& hint) {
    return field(path, [&](forms::FieldSpec& spec) {
        spec.label = label;
        if (spec.hint.isEmpty()) {
            spec.hint = hint;
        }
    });
}

}  // namespace mdbridge::app::pages
