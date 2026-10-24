#pragma once

#include <QTabWidget>

class QFormLayout;
class QLabel;
class QTreeWidget;

namespace mdbridge::app {

struct OpenedProject;

// Sidebar: the project's summary and settings, and its step log.
class ProjectPanel : public QTabWidget {
    Q_OBJECT

public:
    explicit ProjectPanel(QWidget* parent = nullptr);

    void showProject(const OpenedProject& project);

private:
    QLabel* m_placeholder = nullptr;
    QWidget* m_summary = nullptr;
    QFormLayout* m_summaryRows = nullptr;
    QTreeWidget* m_settings = nullptr;
    QTreeWidget* m_log = nullptr;
};

}  // namespace mdbridge::app
