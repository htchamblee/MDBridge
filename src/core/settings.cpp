#include "mdbridge/core/settings.hpp"

#include <algorithm>
#include <cstdlib>
#include <optional>
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

#ifndef _WIN32
// <XDG config home>: $XDG_CONFIG_HOME as given, else ~/.config
fs::path xdgConfigHome() {
    std::string root = environment("XDG_CONFIG_HOME");  // used as given, even relative (as MDATK)
    if (root.empty()) {
        root = posixpath::join(posixpath::expanduser("~"), ".config");
    }
    return pathFromUtf8(root);
}
#endif

}  // namespace

fs::path mdbridgeConfigDir() {
#ifdef _WIN32
    const std::string appData = environment("APPDATA");
    if (!appData.empty()) {
        return pathFromUtf8(appData) / "MDBridge";
    }
    return pathFromUtf8(environment("USERPROFILE")) / "AppData" / "Roaming" / "MDBridge";
#else
    return xdgConfigHome() / "mdbridge";
#endif
}

fs::path userSettingsPath() {
#ifdef _WIN32
    return mdbridgeConfigDir() / "settings.toml";
#else
    return xdgConfigHome() / "mdatk" / "settings.toml";
#endif
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

namespace {

// projectfile._settings_to_change of MDATK 0.1.2 (AP-17): the user settings, read to change one and
// write them all back; nothing when the file is missing. A file that is there but is no TOML (cut
// off, edited by hand, not UTF-8) is never written over, which would lose every setting in it: it is
// kept as settings.toml.unreadable, with a warning, and a new file starts. A file that cannot be read
// at all throws (fs::filesystem_error), so nothing is written.
Table settingsToChange(const fs::path& target) {
    std::string reason;
    try {
        Value data = parseToml(detail::readTextFile(target), detail::displayPath(target));
        if (Table* table = data.get<Table>()) {
            return std::move(*table);
        }
        return {};
    } catch (const fs::filesystem_error& error) {
        if (error.code() == std::errc::no_such_file_or_directory) {
            return {};
        }
        throw;
    } catch (const ProjectError& error) {  // not UTF-8, not TOML
        reason = error.what();
    }
    const std::string shown = detail::displayPath(target);
    if (reason.starts_with(shown + ": ")) {
        reason.erase(0, shown.size() + 2);
    }
    fs::path kept = target;
    kept += ".unreadable";
    fs::rename(target, kept);
    logWarning(shown + " cannot be read (" + reason + "): it is kept as " + detail::displayPath(kept) +
               " and a new settings file starts");
    return {};
}

// projectfile._SETTING_DOCS of MDATK 0.0.14
std::optional<std::string> settingDoc(const std::vector<std::string>& path, const std::vector<const Value*>&) {
    const std::string key = path.empty() ? std::string() : path.back();
    if (key == "recent") {
        return "Recently used project files, newest first.";
    }
    if (key == "cat_proof") {
        return "RNAdo's Cat-proof mode: the keyboard is ignored (untick it with the mouse).";
    }
    return std::nullopt;
}

}  // namespace

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
        Table data = settingsToChange(target);  // the other settings stay (MDATK 0.0.14, 0.1.2)
        const auto entry = std::find_if(data.begin(), data.end(), [](const auto& e) { return e.first == "recent"; });
        if (entry != data.end()) {
            entry->second = std::move(recent);
        } else {
            data.emplace_back("recent", std::move(recent));
        }
        const std::string text = dumpToml(Value(std::move(data)), settingDoc,
                                          {"Molecular Dynamic Analysis Tool Kit (MDATK) user settings."});
        detail::replaceFile(target, text);  // atomic: MDATK writes the file in place
    } catch (const fs::filesystem_error& error) {
        logWarning(std::string("could not update the recent projects list: ") + error.what());
    }
}

}  // namespace mdbridge::core
