#include "render_backend_qt.hpp"

#include <string>

#include <QGuiApplication>
#include <QtGlobal>

#include "mdbridge/render/figure_qt.hpp"
#include "mdbridge/render/text_qt.hpp"
#include "mdbridge/render/video_qt.hpp"
#include "video/encoder.hpp"

namespace mdbridge::app {

RenderBackend qtRenderBackend(const char* program) {
    RenderBackend backend;
    const std::string name = program ? program : "mdbridge";
    backend.start = [name]() -> std::shared_ptr<void> {
        if (QGuiApplication::instance() != nullptr) {
            return {};
        }
        return std::shared_ptr<QGuiApplication>(render::startBatchApplication(name.c_str()).release());
    };
    backend.fonts = []() -> const render::FontSet& { return render::bundledFontSet(); };
    backend.writeFigure = [](const render::Figure& figure, const std::filesystem::path& path, const std::string& format) {
        render::writeFigure(figure, path, format);
    };
    backend.renderVideo = [](const core::Dataset& dataset, const core::Project& project,
                             const std::filesystem::path& path, const core::Reporter& reporter) {
        return render::renderVideo(dataset, project, path, reporter);
    };
    backend.environment = []() {
        return core::Table{{"qt", std::string(qVersion())}, {"ffmpeg", render::video::ffmpegVersionLine()}};
    };
    return backend;
}

}  // namespace mdbridge::app
