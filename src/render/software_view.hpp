#pragma once
// The software view: QPainter drawing for machines without Vulkan 1.3 (older GPUs, virtual
// machines, remote desktops).

#include <QBasicTimer>
#include <QElapsedTimer>
#include <QRasterWindow>

#include "mdbridge/render/view.hpp"

namespace mdbridge::render {

class SoftwareWindow : public QRasterWindow {
    Q_OBJECT

public:
    explicit SoftwareWindow(QString note);

    // The current frame drawn into an image of `size` device-independent pixels.
    QImage render(QSize size, qreal devicePixelRatio) const;

signals:
    void frameRendered();

protected:
    void paintEvent(QPaintEvent* event) override;
    void timerEvent(QTimerEvent* event) override;

private:
    void draw(QPainter& painter, QSizeF size) const;

    QString m_note;
    QElapsedTimer m_clock;
    QBasicTimer m_timer;
};

class SoftwareView : public View {
    Q_OBJECT

public:
    explicit SoftwareView(const QString& reason);
    QImage grabFrame() override;
};

}  // namespace mdbridge::render
