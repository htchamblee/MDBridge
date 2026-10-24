#include "gui.hpp"

#include <iostream>
#include <unordered_set>

#include <QApplication>
#include <QIcon>
#include <QImage>
#include <QTimer>

#include "main_window.hpp"
#include "mdbridge/core/json.hpp"
#include "mdbridge/core/paths.hpp"
#include "mdbridge/version.hpp"

namespace mdbridge::app {

namespace {

constexpr int kScreenshotTimeoutMs = 60000;

std::size_t distinctColours(const QImage& image) {
    const QImage rgb = image.convertToFormat(QImage::Format_RGB32);
    std::unordered_set<QRgb> colours;
    for (int y = 0; y < rgb.height(); ++y) {
        const auto* line = reinterpret_cast<const QRgb*>(rgb.constScanLine(y));
        colours.insert(line, line + rgb.width());
    }
    return colours.size();
}

// Saves the current view's frame and prints what was rendered, as JSON on one line.
int takeScreenshot(MainWindow& window, const GuiOptions& options) {
    render::View* view = window.view();
    const QImage image = view->grabFrame();
    const QString path = QString::fromStdString(options.screenshot);
    const bool saved = !image.isNull() && image.save(path, "PNG");
    const core::Value report = core::Table{
        {"renderer", render::rendererName(view->renderer())},
        {"description", view->description().toStdString()},
        {"fallback_reason", window.fallbackReason().toStdString()},
        {"frames", view->frames()},
        {"width", image.width()},
        {"height", image.height()},
        {"colours", static_cast<std::int64_t>(image.isNull() ? 0 : distinctColours(image))},
        {"image", saved ? core::Value(options.screenshot) : core::Value()},
    };
    std::cout << core::toJson(report, -1) << std::endl;
    if (!saved) {
        std::cerr << "mdbridge: could not " << (image.isNull() ? "read back the 3D view" : "save " + options.screenshot)
                  << '\n';
        return 2;
    }
    return 0;
}

}  // namespace

int runGui(const GuiOptions& options, int& argc, char** argv) {
    QApplication::setApplicationName(QStringLiteral("MDBridge"));
    QApplication::setApplicationVersion(QString::fromLatin1(kVersion));
    QApplication::setDesktopFileName(QStringLiteral("mdbridge"));
    QApplication app(argc, argv);
    QApplication::setWindowIcon(QIcon(QStringLiteral(":/mdbridge/icons/mdbridge.png")));

    const bool checking = !options.screenshot.empty() || options.quitAfter > 0;
    MainWindow window(options.renderer);
    if (!options.project.empty() && !window.openProject(core::pathFromUtf8(options.project), !checking) && checking) {
        return 1;
    }
    window.resize(1280, 800);
    window.show();

    if (!options.screenshot.empty()) {
        // queued: the Vulkan view reports frames from inside its frame callback, where it cannot be read back
        bool taken = false;
        QObject::connect(
            &window, &MainWindow::viewFrameRendered, &app,
            [&](int frames) {
                if (!taken && frames >= options.screenshotFrames) {
                    taken = true;
                    QApplication::exit(takeScreenshot(window, options));
                }
            },
            Qt::QueuedConnection);
        QTimer::singleShot(kScreenshotTimeoutMs, &app, [&] {
            if (!taken) {
                std::cerr << "mdbridge: the 3D view rendered no frame in " << kScreenshotTimeoutMs / 1000 << " s\n";
                QApplication::exit(3);
            }
        });
    }
    if (options.quitAfter > 0) {
        QTimer::singleShot(static_cast<int>(options.quitAfter * 1000), &app, &QApplication::quit);
    }
    return QApplication::exec();
}

}  // namespace mdbridge::app
