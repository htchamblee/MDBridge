#include "cli.hpp"

#include <algorithm>
#include <iostream>
#include <string>

#include <CLI/CLI.hpp>

#ifdef _WIN32
#include <windows.h>
#endif

#include "gui/gui.hpp"
#include "mdbridge/core/json.hpp"
#include "mdbridge/core/paths.hpp"
#include "mdbridge/render/vulkan_probe.hpp"
#include "mdbridge/version.hpp"
#include "project_summary.hpp"

namespace mdbridge::app {

namespace {

const char* kDescription = "MDBridge: analysis and visualisation of MD trajectories of RNA systems (successor of MDATK)";

std::string versionText(const char* program) {
    return std::string(program) + " " + kVersion + " (reads MDATK " + kMdatkReference + " project files)";
}

struct GuiArguments {
    GuiOptions options;
    std::string renderer = "auto";
};

void addGuiOptions(CLI::App& command, GuiArguments& args) {
    command.add_option("project", args.options.project,
                       "project file to open (.mdatk.toml, or a run's MDATK_params_log.toml)");
    command.add_option("--renderer", args.renderer, "3D view: auto (Vulkan 1.3 when available), vulkan or software")
        ->check(CLI::IsMember({"auto", "vulkan", "software"}))
        ->capture_default_str();
    command.add_option("--screenshot", args.options.screenshot,
                       "check: save the 3D view to this PNG after a few frames, print a JSON report and quit");
    command.add_option("--screenshot-frames", args.options.screenshotFrames, "frames to render before --screenshot")
        ->check(CLI::PositiveNumber)
        ->capture_default_str();
    command.add_option("--quit-after", args.options.quitAfter, "check: quit after this many seconds")
        ->check(CLI::NonNegativeNumber);
}

int startGui(GuiArguments& args, int& argc, char** argv) {
    args.options.renderer = *render::rendererFromName(args.renderer);
    return runGui(args.options, argc, argv);  // Qt reads its own options (-platform ...) from argv
}

int projectShow(const std::string& file, bool json) {
    try {
        const OpenedProject project = openProject(core::pathFromUtf8(file));
        if (json) {
            std::cout << core::toJson(projectData(project)) << '\n';
            return 0;
        }
        const std::vector<SummaryRow> rows = summarize(project);
        std::size_t width = 0;
        for (const SummaryRow& row : rows) {
            width = std::max(width, row.label.size());
        }
        for (const SummaryRow& row : rows) {
            std::cout << row.label << std::string(width + 2 - row.label.size(), ' ') << row.value << '\n';
        }
        return 0;
    } catch (const core::ProjectError& error) {
        std::cerr << "mdbridge: " << error.what() << '\n';
        return 1;
    }
}

int graphics(bool json) {
    const render::VulkanReport report = render::probeVulkan();
    const std::string required = std::to_string(render::kVulkanMajor) + "." + std::to_string(render::kVulkanMinor);
    if (json) {
        core::Array devices;
        for (const render::VulkanDevice& d : report.devices) {
            devices.push_back(core::Table{{"name", d.name},
                                          {"type", d.type},
                                          {"api_version", render::vulkanVersionText(d.apiVersion)},
                                          {"meets_minimum", d.meetsMinimum()}});
        }
        const core::Value data = core::Table{
            {"required", required},
            {"loader_found", report.loaderFound},
            {"loader_version", report.loaderFound ? core::Value(render::vulkanVersionText(report.loaderVersion)) : core::Value()},
            {"devices", devices},
            {"renderer", report.usable() ? "vulkan" : "software"},
            {"device", report.usable() ? core::Value(report.devices[static_cast<std::size_t>(report.chosen)].name) : core::Value()},
            {"problem", report.problem},
        };
        std::cout << core::toJson(data) << '\n';
        return 0;
    }
    std::cout << "Vulkan loader  "
              << (report.loaderFound ? "found, Vulkan " + render::vulkanVersionText(report.loaderVersion) : "not found")
              << '\n';
    for (std::size_t i = 0; i < report.devices.size(); ++i) {
        const render::VulkanDevice& d = report.devices[i];
        std::cout << (i ? "               " : "Devices        ") << d.name << " (" << d.type << ", Vulkan "
                  << render::vulkanVersionText(d.apiVersion) << ")"
                  << (static_cast<int>(i) == report.chosen ? "  <- 3D view" : d.meetsMinimum() ? "" : "  needs Vulkan " + required)
                  << '\n';
    }
    if (report.usable()) {
        std::cout << "3D view        Vulkan " << required << " on "
                  << report.devices[static_cast<std::size_t>(report.chosen)].name
                  << " (where the display supports Vulkan; otherwise the software view)\n";
    } else {
        std::cout << "3D view        software view: " << report.problem << '\n';
    }
    return 0;
}

}  // namespace

int cliMain(int argc, char** argv) {
#ifdef _WIN32
    SetConsoleOutputCP(CP_UTF8);
#endif
    CLI::App app(kDescription, "mdbridge");
    char** utf8 = app.ensure_utf8(argv);
    app.set_version_flag("--version", versionText("mdbridge"));
    app.require_subcommand(0, 1);
    app.footer("Without a command the AppImage starts the interface. Licence: all rights reserved (see LICENSE).");

    GuiArguments guiArgs;
    CLI::App* gui = app.add_subcommand("gui", "start the interface");
    addGuiOptions(*gui, guiArgs);

    CLI::App* project = app.add_subcommand("project", "project files (.mdatk.toml)");
    project->require_subcommand(1);
    CLI::App* show = project->add_subcommand("show", "summarise a project file (--json: everything read from it)");
    std::string file;
    bool showJson = false;
    show->add_option("file", file, "project file, or a run's MDATK_params_log.toml")->required();
    show->add_flag("--json", showJson, "print the header, settings and step log as JSON");

    CLI::App* graphicsCommand =
        app.add_subcommand("graphics", "the 3D view's renderer: Vulkan " + std::to_string(render::kVulkanMajor) + "." +
                                           std::to_string(render::kVulkanMinor) + " devices or the software view");
    bool graphicsJson = false;
    graphicsCommand->add_flag("--json", graphicsJson, "print the report as JSON");

    try {
        app.parse(argc, utf8);
    } catch (const CLI::ParseError& error) {
        return app.exit(error);
    }
    if (gui->parsed()) {
        return startGui(guiArgs, argc, argv);
    }
    if (show->parsed()) {
        return projectShow(file, showJson);
    }
    if (graphicsCommand->parsed()) {
        return graphics(graphicsJson);
    }
    std::cout << app.help();
    return 0;
}

int guiEntryMain(int argc, char** argv) {
    CLI::App app(kDescription, "mdbridge-gui");
    char** utf8 = app.ensure_utf8(argv);
    app.set_version_flag("--version", versionText("mdbridge-gui"));
    GuiArguments guiArgs;
    addGuiOptions(app, guiArgs);
    try {
        app.parse(argc, utf8);
    } catch (const CLI::ParseError& error) {
        return app.exit(error);
    }
    return startGui(guiArgs, argc, argv);
}

}  // namespace mdbridge::app
