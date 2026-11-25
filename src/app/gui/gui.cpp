#include "gui.hpp"

#include <cstdlib>
#include <iostream>
#include <unordered_set>

#include <QApplication>
#include <QIcon>
#include <QImage>
#include <QPaintDeviceWindow>
#include <QThreadPool>
#include <QTimer>

#include "main_window.hpp"
#include "mdbridge/core/json.hpp"
#include "mdbridge/core/paths.hpp"
#include "mdbridge/version.hpp"
#include "preview_widget.hpp"

namespace mdbridge::app {

namespace {

constexpr int kScreenshotTimeoutMs = 120000;
constexpr int kRetryMs = 50;  // asks the view for another frame while waiting

std::size_t distinctColours(const QImage& image) {
    const QImage rgb = image.convertToFormat(QImage::Format_RGB32);
    std::unordered_set<QRgb> colours;
    for (int y = 0; y < rgb.height(); ++y) {
        const auto* line = reinterpret_cast<const QRgb*>(rgb.constScanLine(y));
        colours.insert(line, line + rgb.width());
    }
    return colours.size();
}

// The part of the pixels that differ from `background` by more than a little in some channel: how
// much of the view the drawing covers.
double coverage(const QImage& image, QRgb background) {
    const QImage rgb = image.convertToFormat(QImage::Format_RGB32);
    std::size_t covered = 0;
    for (int y = 0; y < rgb.height(); ++y) {
        const auto* line = reinterpret_cast<const QRgb*>(rgb.constScanLine(y));
        for (int x = 0; x < rgb.width(); ++x) {
            const QRgb p = line[x];
            if (std::abs(qRed(p) - qRed(background)) > 8 || std::abs(qGreen(p) - qGreen(background)) > 8 ||
                std::abs(qBlue(p) - qBlue(background)) > 8) {
                ++covered;
            }
        }
    }
    const auto total = static_cast<double>(rgb.width()) * rgb.height();
    return total > 0 ? static_cast<double>(covered) / total : 0.0;
}

core::Value imageReport(const QImage& image, const std::string& path, bool saved) {
    return core::Table{
        {"width", image.width()},
        {"height", image.height()},
        {"colours", static_cast<std::int64_t>(image.isNull() ? 0 : distinctColours(image))},
        {"image", saved ? core::Value(path) : core::Value()},
    };
}

// The 3D view (and the preview) saved, and what was rendered, as JSON on one line.
int takeScreenshot(MainWindow& window, const GuiOptions& options) {
    render::View* view = window.view();
    const QImage image = view->grabFrame();
    const bool saved = !image.isNull() && image.save(QString::fromStdString(options.screenshot), "PNG");
    QRgb background = image.isNull() ? 0 : image.pixel(0, 0);
    if (const auto frame = view->moleculeFrame()) {
        background = qRgb(static_cast<int>(frame->background[0] * 255.0f + 0.5f),
                          static_cast<int>(frame->background[1] * 255.0f + 0.5f),
                          static_cast<int>(frame->background[2] * 255.0f + 0.5f));
    }
    core::Table report{
        {"renderer", render::rendererName(view->renderer())},
        {"description", view->description().toStdString()},
        {"fallback_reason", window.fallbackReason().toStdString()},
        {"frames", view->frames()},
        {"width", image.width()},
        {"height", image.height()},
        {"colours", static_cast<std::int64_t>(image.isNull() ? 0 : distinctColours(image))},
        {"image", saved ? core::Value(options.screenshot) : core::Value()},
        {"molecule", view->moleculeShown()},
        {"coverage", image.isNull() ? 0.0 : coverage(image, background)},
        {"structure", window.structureSummary().toStdString()},
        {"load_problem", window.loadProblem().toStdString()},
    };
    bool previewSaved = true;
    {
        const QImage preview = window.preview()->image();
        bool written = false;
        if (!options.previewScreenshot.empty()) {
            written = !preview.isNull() && preview.save(QString::fromStdString(options.previewScreenshot), "PNG");
            previewSaved = written;
        }
        core::Value entry = imageReport(preview, options.previewScreenshot, written);
        entry.get<core::Table>()->emplace_back("problem", window.preview()->problem().toStdString());
        report.emplace_back("preview", std::move(entry));
    }
    std::cout << core::toJson(core::Value(std::move(report)), -1) << std::endl;
    if (!saved) {
        std::cerr << "mdbridge: could not " << (image.isNull() ? "read back the 3D view" : "save " + options.screenshot)
                  << '\n';
        return 2;
    }
    if (!previewSaved) {
        std::cerr << "mdbridge: could not save the figure preview to " << options.previewScreenshot << '\n';
        return 2;
    }
    return 0;
}

// Whether everything a screenshot should show is there: the project loaded, the view showing its
// structure (when it has one) and the preview drawn (or failed), enough frames drawn.
bool screenshotReady(const MainWindow& window, const GuiOptions& options) {
    const render::View* view = window.view();
    if (view->frames() < options.screenshotFrames) {
        return false;
    }
    if (options.project.empty()) {
        return true;
    }
    if (window.loading()) {
        return false;
    }
    if (window.structureLoaded() && !view->moleculeShown()) {
        return false;
    }
    return window.preview()->settled();
}

}  // namespace

int runGui(const GuiOptions& options, int& argc, char** argv) {
    QApplication::setApplicationName(QStringLiteral("MDBridge"));
    QApplication::setApplicationVersion(QString::fromLatin1(kVersion));
    QApplication::setDesktopFileName(QStringLiteral("mdbridge"));
    QApplication app(argc, argv);
    QApplication::setWindowIcon(QIcon(QStringLiteral(":/mdbridge/icons/mdbridge.png")));

    const bool checking = !options.screenshot.empty() || options.quitAfter > 0;
    int code = 0;
    {
        MainWindow window(options.renderer);
        if (!options.project.empty() && !window.openProject(core::pathFromUtf8(options.project), !checking) &&
            checking) {
            return 1;
        }
        window.resize(1280, 800);
        window.show();

        bool taken = false;
        QTimer retry;
        retry.setSingleShot(true);
        retry.setInterval(kRetryMs);
        if (!options.screenshot.empty()) {
            const auto check = [&] {
                if (taken) {
                    return;
                }
                if (screenshotReady(window, options)) {
                    taken = true;
                    QApplication::exit(takeScreenshot(window, options));
                } else if (!retry.isActive()) {
                    retry.start();
                }
            };
            // a molecule is drawn on demand: ask for frames until everything is there
            QObject::connect(&retry, &QTimer::timeout, &app, [&taken, &window, check] {
                if (taken) {
                    return;
                }
                QWindow* w = window.view()->window();
                if (auto* painted = qobject_cast<QPaintDeviceWindow*>(w)) {
                    painted->update();
                } else if (w) {
                    w->requestUpdate();
                }
                check();
            });
            // queued: the Vulkan view reports frames from inside its frame callback, where it cannot be read back
            QObject::connect(&window, &MainWindow::viewFrameRendered, &app, [check](int) { check(); }, Qt::QueuedConnection);
            QObject::connect(&window, &MainWindow::projectLoaded, &app, check, Qt::QueuedConnection);
            QObject::connect(window.preview(), &PreviewWidget::drawn, &app, check, Qt::QueuedConnection);
            QTimer::singleShot(kScreenshotTimeoutMs, &app, [&] {
                if (!taken) {
                    std::cerr << "mdbridge: the 3D view and the preview were not ready in " << kScreenshotTimeoutMs / 1000
                              << " s (frames " << window.view()->frames() << ", loading " << window.loading()
                              << ", structure shown " << window.view()->moleculeShown() << ")\n";
                    QApplication::exit(3);
                }
            });
        }
        if (options.quitAfter > 0) {
            QTimer::singleShot(static_cast<int>(options.quitAfter * 1000), &app, &QApplication::quit);
        }
        code = QApplication::exec();
        window.cancelLoad();
        // the workers (structure load, figure preview, software renders) end before the window goes
        QThreadPool::globalInstance()->waitForDone();
    }
    return code;
}

}  // namespace mdbridge::app
