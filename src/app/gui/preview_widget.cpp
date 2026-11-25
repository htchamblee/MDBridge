#include "preview_widget.hpp"

#include <algorithm>
#include <cmath>
#include <exception>

#include <QPainter>

#include "mdbridge/render/display_list.hpp"
#include "mdbridge/render/figure_qt.hpp"
#include "mdbridge/render/text_qt.hpp"
#include "project_view.hpp"

namespace mdbridge::app {

namespace {
constexpr int kDebounceMs = 120;  // MDATK's canvas: schedule_render(delay_ms=120)
}

PreviewWidget::PreviewWidget(QWidget* parent) : QWidget(parent) {
    setMinimumSize(160, 90);
    m_message = tr("The figure preview appears here once a project is open.");
    m_debounce.setSingleShot(true);
    m_debounce.setInterval(kDebounceMs);
    connect(&m_debounce, &QTimer::timeout, this, &PreviewWidget::start);
}

void PreviewWidget::setContent(std::shared_ptr<const PreviewContent> content) {
    m_content = std::move(content);
    m_message.clear();
    m_problem.clear();
    m_image = QImage();
    ++m_generation;
    if (m_content) {
        m_message = tr("Building preview ...");
        schedule();
    }
    update();
}

void PreviewWidget::setMessage(const QString& message) {
    m_content.reset();
    m_image = QImage();
    ++m_generation;
    m_message = message;
    update();
}

void PreviewWidget::resizeEvent(QResizeEvent* event) {
    QWidget::resizeEvent(event);
    if (m_content) {
        ++m_generation;
        schedule();
    }
}

void PreviewWidget::schedule() {
    m_debounce.start();
}

void PreviewWidget::start() {
    if (!m_content) {
        return;
    }
    if (m_busy) {
        m_pending = true;  // drawn again when the running drawing ends
        return;
    }
    m_pending = false;
    const qreal dpr = devicePixelRatioF();
    const int width = static_cast<int>(std::lround(this->width() * dpr));
    const int height = static_cast<int>(std::lround(this->height() * dpr));
    if (width < kPreviewMinimumSize || height < kPreviewMinimumSize) {
        return;  // MDATK draws nothing in a space under 50 pixels
    }
    m_busy = true;
    struct Result {
        QImage image;
        QString problem;
    };
    m_background.run(
        [content = m_content, width, height, dpr]() {
            Result r;
            try {
                const render::FontSet& fonts = render::bundledFontSet();
                const render::panels::BuiltFigure built =
                    previewFigure(content->project, content->example, fonts, width, height, content->watermark);
                r.image = render::renderImage(render::buildDisplayList(built.figure, fonts));
                r.image.setDevicePixelRatio(dpr);
            } catch (const std::exception& e) {
                r.problem = tr("Cannot draw this layout: %1").arg(QString::fromUtf8(e.what()));
            }
            return r;
        },
        [this, generation = m_generation](Result r) {
            m_busy = false;
            if (generation == m_generation) {
                m_drawn = generation;
                m_image = r.image;
                m_problem = r.problem;
                m_message = r.problem;
                update();
                emit drawn();
            } else {
                m_pending = true;
            }
            if (m_pending) {
                start();
            }
        });
}

void PreviewWidget::paintEvent(QPaintEvent*) {
    QPainter painter(this);
    painter.fillRect(rect(), palette().window());
    if (!m_image.isNull()) {
        const QSizeF shown = m_image.deviceIndependentSize();
        const QSizeF available = size();
        // a size change shows the last image scaled until the new one is drawn
        const double scale = std::min(available.width() / shown.width(), available.height() / shown.height());
        const QSizeF target = m_drawn == m_generation ? shown : shown * scale;
        const QRectF where(QPointF((available.width() - target.width()) / 2, (available.height() - target.height()) / 2),
                           target);
        painter.setRenderHint(QPainter::SmoothPixmapTransform);
        painter.drawImage(where, m_image);
        return;
    }
    if (!m_message.isEmpty()) {
        painter.setPen(palette().color(QPalette::PlaceholderText));
        painter.drawText(rect().adjusted(12, 12, -12, -12), Qt::AlignCenter | Qt::TextWordWrap, m_message);
    }
}

}  // namespace mdbridge::app
