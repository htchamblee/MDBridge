#include "mdbridge/core/settings.hpp"

#include <cstdlib>
#include <string>
#include <system_error>

#include "mdbridge/core/log.hpp"
#include "mdbridge/core/paths.hpp"
#include "mdbridge/core/project_file.hpp"
#include "mdbridge/core/python.hpp"
#include "mdbridge/core/toml_writer.hpp"
#include "project_internal.hpp"

namespace mdbridge::core {

namespace fs = std::filesystem;

namespace {

std::string environment(const char* name) {
#ifdef _WIN32
    std::size_t size = 0;
    if (_wgetenv_s(&size, nullptr, 0, fs::path(name).c_str()) != 0 || size == 0) {
        return {};
    }
    std::wstring value(size, L'\0');
    if (_wgetenv_s(&size, value.data(), size, fs::path(name).c_str()) != 0) {
        return {};
    }
    value.resize(size ? size - 1 : 0);
    return pathToUtf8(fs::path(value));
#else
    const char* value = std::getenv(name);
    return value ? value : "";
#endif
}

}  // namespace

fs::path userSettingsPath() {
    std::string root = environment("XDG_CONFIG_HOME");  // used as given, even relative
    if (root.empty()) {
#ifdef _WIN32
        root = pathToUtf8(pathFromUtf8(environment("USERPROFILE")) / ".config");  // Python's "~"
#else
        root = posixpath::join(posixpath::expanduser("~"), ".config");
#endif
    }
    return pathFromUtf8(root) / "mdatk" / "settings.toml";
}

std::vector<fs::path> recentProjects() {
    std::vector<fs::path> out;
    Value data;
    try {
        data = parseToml(detail::readTextFile(userSettingsPath()), detail::displayPath(userSettingsPath()));
    } catch (const std::exception&) {  // missing, unreadable or not TOML: no list
        return out;
    }
    const Value* recent = data.find("recent");
    const Array* items = recent ? recent->get<Array>() : nullptr;
    if (!items) {
        return out;
    }
    for (const Value& item : *items) {
        const std::string* text = item.get<std::string>();
        std::error_code ec;
        if (text && fs::is_regular_file(pathFromUtf8(*text), ec)) {
            out.push_back(pathFromUtf8(*text));
            if (out.size() == kRecentMax) {
                break;
            }
        }
    }
    return out;
}

void addRecent(const fs::path& project) {
    const std::string path = detail::absolutePath(project);
    Array recent{Value(path)};
    for (const fs::path& other : recentProjects()) {
        if (pathToUtf8(other) != path && recent.size() < kRecentMax) {
            recent.emplace_back(pathToUtf8(other));
        }
    }
    const fs::path target = userSettingsPath();
    try {
        fs::create_directories(target.parent_path());
        const Value data = Table{{"recent", std::move(recent)}};
        const std::string text = dumpToml(
            data,
            [](const std::vector<std::string>&, const std::vector<const Value*>&) -> std::optional<std::string> {
                return "Recently used project files, newest first.";
            },
            {"Molecular Dynamic Analysis Tool Kit (MDATK) user settings."});
        detail::writeTextFile(target, text);
    } catch (const fs::filesystem_error& error) {
        logWarning(std::string("could not update the recent projects list: ") + error.what());
    }
}

}  // namespace mdbridge::core
