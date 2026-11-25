#include "fonts_command.hpp"

#include <iomanip>
#include <iostream>
#include <sstream>

#include <QFontInfo>
#include <QGuiApplication>
#include <QImage>

#include "mdbridge/core/json.hpp"
#include "mdbridge/render/text_qt.hpp"

namespace mdbridge::app {

namespace {

// strings measured for the report: what MDATK's figures show (labels, tick labels, legends)
constexpr const char* kSamples[] = {"Time (ns)", "RMSD (Å)", "Mg²⁺", "−0.5", "frame 120", "AUGC"};

std::string fixed(double v) {
    std::ostringstream s;
    s << std::fixed << std::setprecision(6) << v;
    return s.str();
}

}  // namespace

int fontsCommand(const char* program, bool json, const std::string& renderPath) {
    try {
        const std::unique_ptr<QGuiApplication> app = render::startBatchApplication(program);
        const render::FontSet& set = render::bundledFontSet();
        const QStringList families = render::registerBundledFonts();
        const QImage image = render::textSampleImage();
        if (!renderPath.empty() && !image.save(QString::fromStdString(renderPath), "PNG")) {
            std::cerr << "mdbridge: cannot write " << renderPath << '\n';
            return 1;
        }
        const std::string pixels = render::imagePixelsSha256(image);
        const QString platform = QGuiApplication::platformName();

        core::Array faces;
        core::Array metrics;
        for (const render::BundledFont& f : render::bundledFonts()) {
            const render::Font& font = set.font(f.face);
            const QFontInfo info(render::figureFont(f.face, 13.888888, image));
            faces.push_back(core::Table{{"file", std::string(f.file)},
                                        {"family", font.family()},
                                        {"style", font.subfamily()},
                                        {"version", font.version()},
                                        {"sha256", std::string(f.sha256)},
                                        {"qt_family", info.family().toStdString()}});
            for (const char* text : kSamples) {
                const render::TextLayout t = set.layout(f.face, text, 10.0, 100.0);
                metrics.push_back(core::Table{{"file", std::string(f.file)},
                                              {"text", std::string(text)},
                                              {"points", 10.0},
                                              {"dpi", 100.0},
                                              {"width", t.width()},
                                              {"height", t.height()},
                                              {"descent", t.descent()}});
            }
        }
        if (json) {
            core::Array qtFamilies;
            for (const QString& name : families) {
                qtFamilies.push_back(name.toStdString());
            }
            const core::Value report = core::Table{
                {"platform", platform.toStdString()},
                {"faces", faces},
                {"qt_families", qtFamilies},
                {"metrics", metrics},
                {"image", core::Table{{"width", image.width()},
                                      {"height", image.height()},
                                      {"pixels_sha256", pixels},
                                      {"file", renderPath.empty() ? core::Value() : core::Value(renderPath)}}},
            };
            std::cout << core::toJson(report) << '\n';
            return 0;
        }
        std::cout << "Figure fonts   bundled DejaVu " << render::kDejaVuVersion << " (Matplotlib's), Qt platform "
                  << platform.toStdString() << '\n';
        for (const render::BundledFont& f : render::bundledFonts()) {
            const render::Font& font = set.font(f.face);
            const render::TextLayout t = set.layout(f.face, kSamples[0], 10.0, 100.0);
            std::cout << "               " << std::left << std::setw(26) << (font.family() + " " + font.subfamily())
                      << std::setw(26) << f.file << "'" << kSamples[0] << "' at 10 pt, 100 dpi: " << fixed(t.width())
                      << " px\n";
        }
        std::cout << "Text image     " << image.width() << " x " << image.height() << ", pixels sha256 " << pixels
                  << (renderPath.empty() ? std::string() : ", written to " + renderPath) << '\n';
        return 0;
    } catch (const render::FontError& error) {
        std::cerr << "mdbridge: " << error.what() << '\n';
        return 1;
    }
}

}  // namespace mdbridge::app
