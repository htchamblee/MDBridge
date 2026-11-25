#pragma once
// Per-user settings (owner decision D16 of milestone 5).
//
// Linux (and other POSIX systems): the recently used project files are in the settings file MDATK
// uses too, $XDG_CONFIG_HOME/mdatk/settings.toml (else ~/.config/mdatk/settings.toml), so both
// programs share the list; settings of MDBridge alone go to its own folder,
// $XDG_CONFIG_HOME/mdbridge (else ~/.config/mdbridge). XDG_CONFIG_HOME is used as given, as MDATK
// does.
//
// Windows: everything is in %APPDATA%\MDBridge (settings.toml holds the recent list there). MDATK's
// Windows build keeps its file in %USERPROFILE%\.config\mdatk, which is not shared.

#include <cstddef>
#include <filesystem>
#include <vector>

namespace mdbridge::core {

inline constexpr std::size_t kRecentMax = 10;

// The folder of MDBridge's own settings (see above; on Windows without APPDATA,
// %USERPROFILE%\AppData\Roaming\MDBridge). Not created here.
std::filesystem::path mdbridgeConfigDir();

// The file that holds the recent list: <XDG config>/mdatk/settings.toml (shared with MDATK);
// on Windows <mdbridgeConfigDir()>/settings.toml.
std::filesystem::path userSettingsPath();

// Recently used projects, newest first: at most kRecentMax that still exist. Never throws.
std::vector<std::filesystem::path> recentProjects();

// Puts `project` (made absolute) first in the list and writes the settings file, keeping every
// other entry of the file as it was (MDATK 0.0.14 keeps RNAdo's cat_proof there; MDATK 0.0.4.1 and
// MDBridge before wrote the list alone). The file is written to <file>.tmp and renamed over the old
// one, so a crash never leaves half a file. A file that cannot be read or is not TOML is replaced
// (as MDATK). A failure is logged as a warning ("could not update the recent projects list:
// <reason>").
void addRecent(const std::filesystem::path& project);

}  // namespace mdbridge::core
