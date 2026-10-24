#include "project_panel.hpp"

#include <QFormLayout>
#include <QHeaderView>
#include <QLabel>
#include <QTreeWidget>
#include <QVBoxLayout>

#include "../project_summary.hpp"

namespace mdbridge::app {

namespace {

using core::Value;

QString toQString(const std::string& text) { return QString::fromUtf8(text.data(), static_cast<qsizetype>(text.size())); }

// One row per entry; tables and arrays of tables become branches.
void addValue(QTreeWidgetItem* parent, const QString& key, const Value& value) {
    auto* item = new QTreeWidgetItem(parent, {key});
    if (const core::Table* table = value.get<core::Table>()) {
        for (const auto& [k, v] : *table) {
            addValue(item, toQString(k), v);
        }
        return;
    }
    const core::Array* array = value.get<core::Array>();
    if (array && !array->empty() && (*array)[0].isTable()) {
        item->setText(1, QObject::tr("%n item(s)", nullptr, static_cast<int>(array->size())));
        for (std::size_t i = 0; i < array->size(); ++i) {
            addValue(item, QString::number(i + 1), (*array)[i]);
        }
        return;
    }
    const QString text = toQString(valueText(value));
    item->setText(1, text.section(QLatin1Char('\n'), 0, 0) + (text.contains(QLatin1Char('\n')) ? QStringLiteral(" ...") : QString()));
    item->setToolTip(1, text);
}

}  // namespace

ProjectPanel::ProjectPanel(QWidget* parent) : QTabWidget(parent) {
    auto* projectTab = new QWidget(this);
    auto* layout = new QVBoxLayout(projectTab);
    m_placeholder = new QLabel(tr("Open a project file (.mdatk.toml) with File > Open project, or drop it on this window."),
                               projectTab);
    m_placeholder->setWordWrap(true);
    m_placeholder->setAlignment(Qt::AlignCenter);
    layout->addWidget(m_placeholder);

    m_summary = new QWidget(projectTab);
    m_summaryRows = new QFormLayout(m_summary);
    m_summaryRows->setContentsMargins(0, 0, 0, 0);
    layout->addWidget(m_summary);

    m_settings = new QTreeWidget(projectTab);
    m_settings->setHeaderLabels({tr("Setting"), tr("Value")});
    m_settings->header()->setSectionResizeMode(0, QHeaderView::ResizeToContents);
    layout->addWidget(m_settings, 1);
    m_summary->hide();
    m_settings->hide();
    addTab(projectTab, tr("Project"));

    m_log = new QTreeWidget(this);
    m_log->setHeaderLabels({tr("Time"), tr("Step"), tr("Title"), tr("Status"), tr("Duration (s)"), tr("Summary")});
    m_log->setRootIsDecorated(false);
    addTab(m_log, tr("Step log"));
}

void ProjectPanel::showProject(const OpenedProject& project) {
    while (m_summaryRows->rowCount() > 0) {
        m_summaryRows->removeRow(0);
    }
    for (const SummaryRow& row : summarize(project)) {
        auto* value = new QLabel(toQString(row.value), m_summary);
        value->setTextFormat(Qt::PlainText);
        value->setWordWrap(true);
        value->setTextInteractionFlags(Qt::TextSelectableByMouse);
        m_summaryRows->addRow(toQString(row.label) + QLatin1Char(':'), value);
    }

    m_settings->clear();
    for (const auto& [key, value] : *project.file.settings.get<core::Table>()) {
        addValue(m_settings->invisibleRootItem(), toQString(key), value);
    }

    m_log->clear();
    for (const Value& entry : project.log) {
        auto field = [&entry](const char* key) {
            const Value* v = entry.find(key);
            return v ? toQString(valueText(*v)) : QString();
        };
        auto* item = new QTreeWidgetItem(m_log, {field("time"), field("step"), field("title"), field("status"),
                                                 field("duration_s"), field("summary").section(QLatin1Char('\n'), 0, 0)});
        if (field("status") == QLatin1String("failed")) {
            item->setForeground(3, QColor(0xc6, 0x28, 0x28));
        }
        const QString output = field("output");
        item->setToolTip(5, field("summary") + (output.isEmpty() ? QString() : QStringLiteral("\n\n") + output.right(2000)));
    }
    for (int column = 0; column < m_log->columnCount(); ++column) {
        m_log->resizeColumnToContents(column);
    }
    setTabText(1, project.log.empty() ? tr("Step log") : tr("Step log (%1)").arg(project.log.size()));

    m_placeholder->hide();
    m_summary->show();
    m_settings->show();
    setCurrentIndex(0);
}

}  // namespace mdbridge::app
