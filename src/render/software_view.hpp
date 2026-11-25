#pragma once
// The software view: QPainter drawing for machines without Vulkan 1.3 (older GPUs, virtual
// machines, remote desktops). A molecule is drawn by the built-in software renderer
// (render/molecule/software.hpp, the structure panel's renderer) on a worker thread: the window
// keeps showing the last image (scaled while the window is resized) until the next one is ready,
// and a request made while one renders replaces any still waiting, so dragging never queues work.

#include <cstdint>
#include <memory>
#include <mutex>

#include <QBasicTimer>
#include <QElapsedTimer>
#include <QImage>
#include <QRasterWindow>

#include "mdbridge/render/view.hpp"

namespace mdbridge::render {

class SoftwareWindow : public QRasterWindow {
    Q_OBJECT

public:
    explicit SoftwareWindow(QString note);
    ~SoftwareWindow() override;

    // The current frame drawn into an image of `size` device-independent pixels (the molecule as
    // last rendered, scaled to `size`).
    QImage render(QSize size, qreal devicePixelRatio) const;

    // Draws `frame` through `camera` from now on (nullptr: the placeholder).
    void setMolecule(std::shared_ptr<const MoleculeFrame> frame, const ViewCamera& camera);
    // The image shown is the current molecule and camera at the window's size.
    bool moleculeShown() const;

signals:
    void frameRendered();

protected:
    void paintEvent(QPaintEvent* event) override;
    void timerEvent(QTimerEvent* event) override;
    void resizeEvent(QResizeEvent* event) override;

private:
    struct Receiver;  // the worker's way back to this window (cleared when it is destroyed)

    void draw(QPainter& painter, QSizeF size) const;
    void drawPlaceholder(QPainter& painter, QSizeF size) const;
    void scheduleRender();
    void startRender();
    void rendered(std::uint64_t generation, const QImage& image, const QString& problem);

    QString m_note;
    QElapsedTimer m_clock;
    QBasicTimer m_timer;

    std::shared_ptr<const MoleculeFrame> m_frame;
    ViewCamera m_camera;
    std::uint64_t m_generation = 0;  // of the molecule, camera and size requested
    std::uint64_t m_imageGeneration = 0;
    QImage m_image;    // the last rendered molecule (device pixels)
    QString m_problem;  // why the last render failed
    bool m_busy = false;
    bool m_pending = false;
    std::shared_ptr<Receiver> m_receiver;
};

class SoftwareView : public View {
    Q_OBJECT

public:
    explicit SoftwareView(const QString& reason);
    QImage grabFrame() override;
    bool moleculeShown() const override;

protected:
    bool showMolecule(std::shared_ptr<const MoleculeFrame> frame, const ViewCamera& camera) override;
};

}  // namespace mdbridge::render
