#pragma once
// The shared libraries this process has loaded, as the dynamic loader opened them (full paths),
// in load order; empty on platforms where MDBridge does not list them (Windows). `graphics --json`
// and the report of `gui --script` give them, so the smoke test can check that the AppImage's
// fallback copies are used only for the libraries a host lacks (packaging/linux/AppRun).

#include <string>
#include <vector>

namespace mdbridge::app {

std::vector<std::string> loadedLibraries();

}  // namespace mdbridge::app
