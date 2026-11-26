#pragma once
// The read-only figure preview of the main window (MDATK's GUI canvas without its editing): the
// project's figure built with the example dataset and watermarked (app/commands/project_view.hpp
// previewFigure), drawn on a worker thread at the widget's size whenever the content or the size
// changes. The interactive layout editor is milestone 5.

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

#include <QImage>
#include <QTimer>
#include <QWidget>

#include "background.hpp"
#include "mdbridge/core/dataset.hpp"
#include "mdbridge/core/project.hpp"

namespace mdbridge::app {

// What the preview draws: the project, its example dataset and the watermark.
struct PreviewContent {
    core::Project project;
    core::Dataset example{{}, {}};
    std::string watermark;
};

class PreviewWidget : public QWidget {
    Q_OBJECT

public:
    explicit PreviewWidget(QWidget* parent = nullptr);

    // Draws `content` (nullptr: the hint to open a project).
    void setContent(std::shared_ptr<const PreviewContent> content);
    // A message instead of a figure (e.g. "Cannot draw this layout: ...").
    void setMessage(const QString& message);

    // The figure as last drawn (null before the first one), and why the last drawing failed.
    QImage image() const { return m_image; }
    QString problem() const { return m_problem; }
    // The warnings the last drawing logged ("WARNING: ..."; layout warnings, centre-on misses).
    const std::vector<std::string>& warnings() const { return m_warnings; }
    // The image shows the current content at the current size.
    bool current() const { return m_content && !m_image.isNull() && m_drawn == m_generation; }
    // Nothing more to draw: no content (a message), or the current content drawn or failed.
    bool settled() const { return !m_content || (m_drawn == m_generation && !m_busy); }

    QSize sizeHint() const override { return {640, 360}; }

signals:
    // A drawing finished (successfully or not: problem()).
    void drawn();

protected:
    void paintEvent(QPaintEvent* event) override;
    void resizeEvent(QResizeEvent* event) override;

private:
    void schedule();
    void start();

    std::shared_ptr<const PreviewContent> m_content;
    QImage m_image;
    QString m_problem;
    std::vector<std::string> m_warnings;
    QString m_message;
    std::uint64_t m_generation = 0;
    std::uint64_t m_drawn = 0;
    bool m_busy = false;
    bool m_pending = false;
    QTimer m_debounce;
    Background m_background{this};
};

}  // namespace mdbridge::app
