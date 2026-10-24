#pragma once
// User settings shared with MDATK: $XDG_CONFIG_HOME/mdatk/settings.toml (else
// ~/.config/mdatk/settings.toml), which holds the recently used project files.

#include <cstddef>
#include <filesystem>
#include <vector>

namespace mdbridge::core {

inline constexpr std::size_t kRecentMax = 10;

std::filesystem::path userSettingsPath();

// Recently used projects, newest first: at most kRecentMax that still exist. Never throws.
std::vector<std::filesystem::path> recentProjects();

// Puts `project` (made absolute) first in the list and writes the settings file; a failure is
// logged as a warning ("could not update the recent projects list: <reason>").
void addRecent(const std::filesystem::path& project);

}  // namespace mdbridge::core
